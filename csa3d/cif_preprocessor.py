"""mmCIF file preprocessing and validation utilities.

Handles preprocessing of mmCIF (Macromolecular Crystallographic Information File)
files downloaded from PDBe to ensure they can be reliably parsed by BioPython.
Specific issues handled include missing atomic properties (occupancies, B-factors)
and invalid numeric values that can cause parsing failures or incorrect
coordinate interpretation.

This module is automatically used during PdbSite.from_list() when the CIF file
needs validation before parsing.
"""

import tempfile
import re
import os
import logging

logger = logging.getLogger(__name__)


def cleanup_temp_cif(cif_path: str, original_path: str) -> None:
    """Remove temporary CIF file if one was created during validation.
    
    Safe cleanup of preprocessing artifacts. Only removes the file if it's
    a temporary file (differs from the original), preventing accidental
    deletion of original data files.
    
    Args:
        cif_path: Path to the file to potentially clean up. May be a temporary
            file created by validate_cif().
        original_path: Path to the original input file. If cif_path differs from
            this, it's assumed to be temporary and will be removed.
    """
    if cif_path != original_path and os.path.exists(cif_path):
        try:
            os.unlink(cif_path)
            logger.debug(f'Cleaned up temporary CIF: {cif_path}')
        except Exception as e:
            logger.warning(f'Failed to clean up temporary CIF {cif_path}: {e}')


def validate_cif(cif_path: str, occupancy_default: float = 1.0, b_factor_default: float = 20.0) -> str:
    """Validate and fix mmCIF file to ensure BioPython compatibility.
    
    Preprocesses mmCIF files to handle common issues that prevent successful
    parsing or lead to incorrect coordinate interpretation:
    
    - Missing occupancy column: Added with default value (typically 1.0 for
      fully occupied atoms)
    - Missing B-factor column: Added with default value (typically 20.0 Ų)
    - Invalid occupancy values: Fixed to valid range [0, 1]
    - Negative B-factors: Replaced with default (sometimes appears in low-quality
      structures)
    - Missing atom records: Filled with placeholder values
    
    If no issues are detected, returns the original file path unchanged. If
    fixes are needed, writes a corrected version to a temporary file and
    returns that path. Caller should use cleanup_temp_cif() after parsing.
    
    Args:
        cif_path: Path to input mmCIF file from PDBe or other source.
        occupancy_default: Default occupancy value to use for missing or invalid
            occupancy data. Typically 1.0 for fully occupied atoms. Defaults to 1.0.
        b_factor_default: Default B-factor (temperature factor) value in Ų for
            missing or invalid data. Typical values range 10-30 depending on
            structure quality. Defaults to 20.0.
    
    Returns:
        str: Path to validated CIF file. Either:
            - The original cif_path if no issues were found
            - Path to a new temporary file if fixes were applied
        
        Caller should call cleanup_temp_cif(result, original_path) after use
        to remove temporary files.
    
    Example:
        >>> from Bio.PDB import MMCIFParser
        >>> cif_path = "original.cif"
        >>> fixed_path = validate_cif(cif_path)
        >>> parser = MMCIFParser()
        >>> structure = parser.get_structure('', fixed_path)
        >>> cleanup_temp_cif(fixed_path, cif_path)
    """
    if not os.path.exists(cif_path):
        logger.warning(f'CIF file not found: {cif_path}')
        return cif_path

    with open(cif_path, 'r') as f:
        lines = f.readlines()

    # Find atom_site loop block
    in_atom_loop = False
    atom_loop_start = -1
    field_names = []
    occupancy_col = -1
    b_factor_col = -1

    for i, line in enumerate(lines):
        if 'loop_' in line and i + 1 < len(lines) and '_atom_site.' in lines[i + 1]:
            in_atom_loop = True
            atom_loop_start = i
            break

    if not in_atom_loop:
        return cif_path

    # Collect field names and find occupancy/B-factor columns
    i = atom_loop_start + 1
    while i < len(lines):
        line = lines[i].strip()
        if line.startswith('_atom_site.'):
            field_name = line.replace('_atom_site.', '').strip()
            field_names.append(field_name)
            if field_name == 'occupancy':
                occupancy_col = len(field_names) - 1
            elif field_name == 'B_iso_or_equiv':
                b_factor_col = len(field_names) - 1
            i += 1
        else:
            break

    data_start = i
    needs_fix = False

    # Copy header and field definitions
    fixed_lines = lines[:data_start]

    # Check if we need to add missing columns
    if occupancy_col == -1:
        logger.info(f'Adding missing occupancy column to {cif_path}')
        field_names.append('occupancy')
        occupancy_col = len(field_names) - 1
        needs_fix = True
        fixed_lines.append('_atom_site.occupancy\n')

    if b_factor_col == -1:
        logger.info(f'Adding missing B_iso_or_equiv column to {cif_path}')
        field_names.append('B_iso_or_equiv')
        b_factor_col = len(field_names) - 1
        needs_fix = True
        fixed_lines.append('_atom_site.B_iso_or_equiv\n')

    # Process atom records - check for invalid values
    replacements = {'occupancy': 0, 'b_factor': 0, 'negative_b_factor': 0}

    for i in range(data_start, len(lines)):
        line = lines[i]
        if not line.strip() or line.strip().startswith('_') or line.strip().startswith('loop_'):
            fixed_lines.extend(lines[i:])
            break

        if not line.strip():
            fixed_lines.append(line)
            continue

        parts = re.findall(r"'[^']*'|\"[^\"]*\"|[^\s]+", line)

        # Extend with missing columns if needed
        while len(parts) < len(field_names):
            if len(parts) == occupancy_col:
                parts.append(str(occupancy_default))
                replacements['occupancy'] += 1
                needs_fix = True
            elif len(parts) == b_factor_col:
                parts.append(str(b_factor_default))
                replacements['b_factor'] += 1
                needs_fix = True
            else:
                parts.append('?')

        # Fix occupancy - check for invalid values
        if occupancy_col >= 0 and occupancy_col < len(parts):
            try:
                occ_val = float(parts[occupancy_col])
                if occ_val < 0 or occ_val > 1.0:
                    logger.debug(f'Fixing invalid occupancy {occ_val}')
                    parts[occupancy_col] = str(occupancy_default)
                    replacements['occupancy'] += 1
                    needs_fix = True
            except ValueError:
                if parts[occupancy_col] in ('?', '.'):
                    parts[occupancy_col] = str(occupancy_default)
                    replacements['occupancy'] += 1
                    needs_fix = True

        # Fix B-factor - check for invalid/negative values
        if b_factor_col >= 0 and b_factor_col < len(parts):
            try:
                b_val = float(parts[b_factor_col])
                if b_val < 0:
                    logger.debug(f'Fixing negative B-factor {b_val}')
                    parts[b_factor_col] = str(b_factor_default)
                    replacements['negative_b_factor'] += 1
                    needs_fix = True
            except ValueError:
                if parts[b_factor_col] in ('?', '.'):
                    parts[b_factor_col] = str(b_factor_default)
                    replacements['b_factor'] += 1
                    needs_fix = True

        fixed_lines.append(' '.join(parts) + '\n')

    if not needs_fix:
        logger.debug(f'No issues found in {cif_path}')
        return cif_path

    logger.info(f'Validated CIF {cif_path}: occupancy={replacements["occupancy"]}, '
               f'b_factor={replacements["b_factor"]}, negative_b_factor={replacements["negative_b_factor"]}')

    # Write to temp file
    tmp = tempfile.NamedTemporaryFile(mode='w', suffix='.cif', delete=False)
    tmp_path = tmp.name
    tmp.writelines(fixed_lines)
    tmp.close()

    return tmp_path
