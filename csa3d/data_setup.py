"""Automatic data setup and download for CSA-3D.

Downloads and organizes required data files from official sources.
Includes bundled datasets that have been retired from the public domain.
"""

import csv
import json
import logging
import os
import shutil
import subprocess
import tarfile
import tempfile
from collections import defaultdict
from pathlib import Path
from typing import Optional
import urllib.request
import urllib.error
from urllib.parse import urljoin

logger = logging.getLogger(__name__)


class DataDownloader:
    """Download and setup CSA-3D data files."""

    # Data source URLs
    SIFTS_URL = "https://www.ebi.ac.uk/pdbe/download/sifts/csv/"
    KEGG_URL = "https://rest.kegg.jp/get/"
    RHEA_URL = "https://www.rhea-db.org/download/rhea/"
    RHEA_FTP_URL = "https://ftp.expasy.org/databases/rhea/ctfiles/"
    RHEA_TSV_URL = "https://ftp.expasy.org/databases/rhea/tsv/"
    PDBE_URL = "https://www.ebi.ac.uk/pdbe/download/"
    MCSA_URL = "https://www.ebi.ac.uk/thornton-srv/m-csa/download/"
    MCSA_API_URL = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/"

    def __init__(self, data_dir: Path):
        """Initialize data downloader.

        Args:
            data_dir: Root directory for CSA-3D data
        """
        self.data_dir = Path(data_dir)
        self.data_dir.mkdir(parents=True, exist_ok=True)
        # Path to bundled data files in the repository
        self.bundled_data_dir = Path(__file__).parent.parent / "data"

    def download_file(self, url: str, output_path: Path) -> bool:
        """Download a file from URL.

        Args:
            url: URL to download from
            output_path: Path to save file

        Returns:
            True if successful, False otherwise
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)

        try:
            logger.info(f"Downloading {url} → {output_path}")
            urllib.request.urlretrieve(url, output_path)
            logger.info(f"✓ Downloaded {output_path.name}")
            return True
        except urllib.error.URLError as e:
            logger.warning(f"✗ Failed to download {url}: {e}")
            return False
        except Exception as e:
            logger.error(f"Unexpected error downloading {url}: {e}")
            return False

    def copy_bundled_file(self, filename: str, output_path: Path) -> bool:
        """Copy a bundled data file from the repository.

        Args:
            filename: Name of file in bundled data directory
            output_path: Path to save file

        Returns:
            True if successful, False otherwise
        """
        output_path.parent.mkdir(parents=True, exist_ok=True)
        source = self.bundled_data_dir / filename

        if not source.exists():
            logger.warning(f"✗ Bundled file not found: {source}")
            return False

        try:
            logger.info(f"Copying {filename} → {output_path}")
            shutil.copy2(source, output_path)
            logger.info(f"✓ Copied {output_path.name}")
            return True
        except Exception as e:
            logger.error(f"Error copying {filename}: {e}")
            return False

    def compile_rhea_reactions(self, temp_dir: Path) -> Optional[Path]:
        """Compile RHEA reactions from downloaded data.

        Downloads RHEA .rd files, EC mappings, and directions, then compiles
        the reaction_components.csv file.

        Args:
            temp_dir: Temporary directory for intermediate files

        Returns:
            Path to compiled CSV file if successful, None otherwise
        """
        try:
            logger.info("Compiling RHEA reaction data from source files...")

            # Download files
            rhea_dir = temp_dir / "rhea_download"
            rhea_dir.mkdir(parents=True, exist_ok=True)

            # Download rhea-rd.tar.gz
            rd_tar = rhea_dir / "rhea-rd.tar.gz"
            if not self.download_file(
                urljoin(self.RHEA_FTP_URL, "rhea-rd.tar.gz"), rd_tar
            ):
                return None

            # Download rhea2ec.tsv
            rhea2ec_file = rhea_dir / "rhea2ec.tsv"
            if not self.download_file(
                urljoin(self.RHEA_TSV_URL, "rhea2ec.tsv"), rhea2ec_file
            ):
                return None

            # Download rhea-directions.tsv
            directions_file = rhea_dir / "rhea-directions.tsv"
            if not self.download_file(
                urljoin(self.RHEA_TSV_URL, "rhea-directions.tsv"), directions_file
            ):
                return None

            # Extract tar file
            logger.info("Extracting RHEA reaction files...")
            with tarfile.open(rd_tar, "r:gz") as tar:
                tar.extractall(rhea_dir)

            # Compile reactions using the compilation script logic
            logger.info("Compiling reaction components...")
            output_csv = rhea_dir / "rhea_reaction_components.csv"

            # Map directional reaction IDs to master IDs
            directional2master = {}
            with open(directions_file, "r") as f:
                next(f)
                for line in csv.reader(f, delimiter="\t"):
                    master = line[0]
                    lr = line[1]
                    rl = line[2]
                    directional2master[lr] = master
                    directional2master[rl] = master

            # Map ECs to master IDs
            master2ec = defaultdict(set)
            with open(rhea2ec_file, "r") as f:
                next(f)
                for line in csv.reader(f, delimiter="\t"):
                    master = line[2]
                    ec = line[3]
                    master2ec[master].add(ec)

            # Extract reactions from .rd files
            reactions = defaultdict(set)
            rd_files_dir = rhea_dir / "rd"
            if rd_files_dir.exists():
                for rd_file in sorted(rd_files_dir.glob("*.rd")):
                    rhea_id = rd_file.stem
                    with open(rd_file, "r") as f:
                        for line in f:
                            if line.startswith("CHEBI:"):
                                chebi_id = line.strip().split(":")[-1]
                                reactions[rhea_id].add(chebi_id)

            # Write output CSV
            with open(output_csv, "w") as f:
                f.write("RHEA_MASTER_ID,RHEA_DIRECTIONAL_ID,EC,CHEBI_IDS\n")
                for k, v in sorted(reactions.items()):
                    master = directional2master.get(k, int(k) - 1)
                    ecs = master2ec.get(str(master), set())
                    for ec in sorted(ecs):
                        f.write(f"{master},{k},{ec},{';'.join(sorted(v))}\n")

            logger.info(f"✓ Compiled {len(reactions)} RHEA reactions")
            return output_csv

        except Exception as e:
            logger.error(f"Error compiling RHEA data: {e}")
            return None

    def compile_kegg_reactions(self) -> Optional[Path]:
        """Compile KEGG reaction components using Biopython REST API.

        Returns:
            Path to compiled CSV file if successful, None otherwise
        """
        try:
            logger.info("Compiling KEGG reaction data from REST API...")

            try:
                from Bio.KEGG import REST
            except ImportError:
                logger.error("Biopython not installed. Cannot compile KEGG data.")
                logger.info("Install with: pip install biopython")
                return None

            # Get EC-compound mappings from KEGG
            kegg_dir = self.data_dir / "kegg"
            kegg_dir.mkdir(parents=True, exist_ok=True)
            output_csv = kegg_dir / "kegg_reaction_components.csv"

            logger.info("Querying KEGG REST API for enzyme-compound mappings...")
            request = REST.kegg_link("compound", "enzyme").readlines()

            ec_compounds = defaultdict(list)
            for line in request:
                line = line.strip().split("\t")
                ec = line[0].split(":")[-1]
                compound = line[1].split(":")[-1]
                ec_compounds[ec].append(compound)

            # Write output CSV
            with open(output_csv, "w") as f:
                f.write("EC,KEGG_COMPOUNDS\n")
                for ec in sorted(ec_compounds.keys()):
                    compounds = sorted(set(ec_compounds[ec]))
                    f.write(f"{ec},{';'.join(compounds)}\n")

            logger.info(f"✓ Compiled {len(ec_compounds)} KEGG EC entries")
            return output_csv

        except Exception as e:
            logger.error(f"Error compiling KEGG data: {e}")
            return None

    def download_mcsa_catalytic_residues(self) -> bool:
        """Download M-CSA catalytic residue data from API.

        Downloads homologues and reference residue annotations from the M-CSA API.

        Returns:
            True if successful
        """
        try:
            logger.info("Downloading M-CSA catalytic residue data from API...")
            mcsa_dir = self.data_dir / "mcsa"
            mcsa_dir.mkdir(parents=True, exist_ok=True)

            # Download homologues residues
            homologues_url = urljoin(self.MCSA_API_URL, "homologues_residues.json")
            homologues_file = mcsa_dir / "catalytic_residues_homologues.json"

            logger.info("Downloading homologues residues...")
            if not self.download_file(homologues_url, homologues_file):
                logger.warning("Failed to download homologues residues")
                return False

            # Download reference residues
            reference_url = urljoin(self.MCSA_API_URL, "references_residues.json")
            reference_file = mcsa_dir / "catalytic_residues_references.json"

            logger.info("Downloading reference residues...")
            if not self.download_file(reference_url, reference_file):
                logger.warning("Failed to download reference residues")
                return False

            logger.info("✓ Downloaded M-CSA catalytic residue data")
            return True

        except Exception as e:
            logger.error(f"Error downloading M-CSA data: {e}")
            return False

    def setup_sifts(self) -> bool:
        """Setup SIFTS (PDB-UniProt mappings).

        Downloads from EBI PDBe SIFTS data.

        Returns:
            True if successful
        """
        logger.info("Setting up SIFTS data...")
        sifts_dir = self.data_dir / "sifts"
        sifts_dir.mkdir(parents=True, exist_ok=True)

        files = {
            "uniprot_pdb.csv": "uniprot_pdb.csv",
            "pdb_chain_uniprot.csv": "pdb_chain_uniprot.csv",
            "pdb_chain_enzyme.csv": "pdb_chain_enzyme.csv",
        }

        success = True
        for filename, remote_name in files.items():
            url = urljoin(self.SIFTS_URL, remote_name)
            output = sifts_dir / filename
            if not self.download_file(url, output):
                success = False

        return success

    def setup_kegg(self) -> bool:
        """Setup KEGG reaction data.

        Attempts to compile fresh KEGG data using Biopython REST API.
        Falls back to bundled data if compilation fails.

        Returns:
            True if successful
        """
        logger.info("Setting up KEGG data...")
        kegg_dir = self.data_dir / "kegg"
        kegg_dir.mkdir(parents=True, exist_ok=True)
        output = kegg_dir / "kegg_reaction_components.csv"

        # Try to compile fresh data
        try:
            compiled_file = self.compile_kegg_reactions()
            if compiled_file and compiled_file.exists():
                return True
        except Exception as e:
            logger.warning(f"Failed to compile fresh KEGG data: {e}")

        # Fall back to bundled data
        logger.info("Falling back to bundled KEGG data...")
        return self.copy_bundled_file("kegg_reaction_components.csv", output)

    def setup_rhea(self) -> bool:
        """Setup RHEA reaction data.

        Attempts to download and compile the latest RHEA data from source.
        Falls back to bundled pre-compiled data if download fails.

        Returns:
            True if successful
        """
        logger.info("Setting up RHEA data...")
        rhea_dir = self.data_dir / "rhea"
        rhea_dir.mkdir(parents=True, exist_ok=True)
        output = rhea_dir / "rhea_reaction_components.csv"

        # Try to download and compile fresh data
        try:
            with tempfile.TemporaryDirectory() as temp_dir:
                compiled_file = self.compile_rhea_reactions(Path(temp_dir))
                if compiled_file and compiled_file.exists():
                    shutil.copy2(compiled_file, output)
                    logger.info(f"✓ Successfully compiled fresh RHEA data")
                    return True
        except Exception as e:
            logger.warning(f"Failed to compile fresh RHEA data: {e}")

        # Fall back to bundled data
        logger.info("Falling back to bundled RHEA data...")
        return self.copy_bundled_file("rhea_reaction_components.csv", output)

    def setup_pdbe(self) -> bool:
        """Setup PDBe data (ions, cofactors, resolution).

        Returns:
            True if successful
        """
        logger.info("Setting up PDBe data...")
        pdbe_dir = self.data_dir / "pdbe"
        pdbe_dir.mkdir(parents=True, exist_ok=True)

        files = {
            "pdb_ions.csv": "pdb_ions.csv",
            "pdb_cofactors.csv": "pdb_cofactors.csv",
            "resolu.idx": "resolu.idx",
        }

        success = True
        for filename, remote_name in files.items():
            url = urljoin(self.PDBE_URL, remote_name)
            output = pdbe_dir / filename
            if not self.download_file(url, output):
                success = False

        return success

    def setup_mcsa(self) -> bool:
        """Setup M-CSA catalytic site data.

        Downloads latest catalytic residue annotations from M-CSA API.
        Uses bundled cofactor data that has been retired from public domain.

        Returns:
            True if successful
        """
        logger.info("Setting up M-CSA data...")
        mcsa_dir = self.data_dir / "mcsa"
        mcsa_dir.mkdir(parents=True, exist_ok=True)

        success = True

        # Download latest catalytic residue data from API
        if not self.download_mcsa_catalytic_residues():
            logger.warning("Failed to download M-CSA catalytic residue data from API")
            success = False

        # Use bundled cofactor data
        output = mcsa_dir / "mcsa_cofactors.csv"
        if not self.copy_bundled_file("mcsa_cofactors.csv", output):
            success = False

        return success

    def setup_metals(self) -> bool:
        """Setup metals reference data.

        Uses bundled metals list that has been retired from the public domain.

        Returns:
            True if successful
        """
        logger.info("Setting up metals reference data...")
        pdbe_dir = self.data_dir / "pdbe"
        pdbe_dir.mkdir(parents=True, exist_ok=True)

        output = pdbe_dir / "metals.csv"
        return self.copy_bundled_file("metals.csv", output)

    def setup_all(self, skip_optional: bool = False) -> dict:
        """Setup all data files.

        Args:
            skip_optional: Skip optional data sources (KEGG)

        Returns:
            Dictionary with setup status for each source
        """
        logger.info(f"Starting automatic CSA-3D data setup to {self.data_dir}")

        results = {}
        results["sifts"] = self.setup_sifts()
        results["rhea"] = self.setup_rhea()
        results["pdbe"] = self.setup_pdbe()
        results["mcsa"] = self.setup_mcsa()
        results["metals"] = self.setup_metals()

        if not skip_optional:
            results["kegg"] = self.setup_kegg()

        # Summary
        logger.info("\n" + "=" * 60)
        logger.info("Data Setup Summary:")
        for source, success in results.items():
            status = "✓ Complete" if success else "✗ Failed/Skipped"
            logger.info(f"  {source.upper()}: {status}")
        logger.info("=" * 60)

        return results


def setup_data(data_dir: Optional[str] = None) -> bool:
    """Setup CSA-3D data automatically.

    Args:
        data_dir: Directory for data files. Defaults to ~/.csa3d/data

    Returns:
        True if setup successful
    """
    if data_dir is None:
        data_dir = str(Path.home() / ".csa3d" / "data")

    logger.info(f"Setting up CSA-3D data in: {data_dir}")

    downloader = DataDownloader(Path(data_dir))
    results = downloader.setup_all()

    success = all(results.values())
    if success:
        logger.info(f"\n✓ Data setup complete!")
        logger.info(f"Set environment variable: export CSA3D_DATA_DIR={data_dir}")
    else:
        logger.warning("\n⚠ Some data sources failed to download")
        logger.info("Please check logs above and retry")

    return success


if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(levelname)s: %(message)s"
    )
    setup_data()
