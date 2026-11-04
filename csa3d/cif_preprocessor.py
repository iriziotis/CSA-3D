"""mmCIF file preprocessing utilities.

Handles preprocessing of mmCIF files to ensure completeness by filling
missing atomic data (occupancies, B-factors) with sensible placeholder values.
"""

import tempfile
import re
import os
import logging

logger = logging.getLogger(__name__)


def cleanup_temp_cif(cif_path: str, original_path: str) -> None:
    """Clean up temporary CIF file if it differs from original.

    Args:
        cif_path: Path to potentially temporary file
        original_path: Path to original file
    """
    if cif_path != original_path and os.path.exists(cif_path):
        try:
            os.unlink(cif_path)
            logger.debug(f'Cleaned up temporary CIF: {cif_path}')
        except Exception as e:
            logger.warning(f'Failed to clean up temporary CIF {cif_path}: {e}')


def validate_cif(cif_path: str, occupancy_default: float = 1.0, b_factor_default: float = 20.0) -> str:
    """Validate and fix CIF file to ensure it can be parsed by BioPython.

    Handles:
    - Missing occupancy columns
    - Missing B-factor columns
    - Invalid numeric values
    - Missing required fields

    Args:
        cif_path: Path to input mmCIF file
        occupancy_default: Default occupancy value
        b_factor_default: Default B-factor value

    Returns:
        Path to validated CIF file (temp if modified, original if valid)
    """
    if not os.path.exists(cif_path):
        logger.warning(f'CIF file not found: {cif_path}')
        return cif_path

    with open(cif_path, 'r') as f:
        lines = f.readlines()

    # Find atom_site loop
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

    # Collect field names
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

    # Always copy everything from beginning up to and including field definitions
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

    # Process atom records - ALWAYS check for invalid values
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
