"""Automatic data setup and download for CSA-3D.

Downloads and organizes required data files from official sources.
Includes bundled datasets that have been retired from the public domain.
"""

import logging
import os
import shutil
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
    PDBE_URL = "https://www.ebi.ac.uk/pdbe/download/"
    MCSA_URL = "https://www.ebi.ac.uk/thornton-srv/m-csa/download/"

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

        Returns:
            True if successful
        """
        logger.info("Setting up KEGG data...")
        kegg_dir = self.data_dir / "kegg"
        kegg_dir.mkdir(parents=True, exist_ok=True)

        # KEGG reactions endpoint
        output = kegg_dir / "kegg_reaction_components.csv"
        logger.info("Note: KEGG data requires special API handling")
        logger.info("Please obtain KEGG reaction data from: https://www.kegg.jp/")
        return False

    def setup_rhea(self) -> bool:
        """Setup RHEA reaction data.

        Uses bundled compiled data that requires the raw RHEA .rd files.

        Returns:
            True if successful
        """
        logger.info("Setting up RHEA data...")
        rhea_dir = self.data_dir / "rhea"
        rhea_dir.mkdir(parents=True, exist_ok=True)

        output = rhea_dir / "rhea_reaction_components.csv"
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

        Uses bundled data files that have been retired from the public domain.

        Returns:
            True if successful
        """
        logger.info("Setting up M-CSA data...")
        mcsa_dir = self.data_dir / "mcsa"
        mcsa_dir.mkdir(parents=True, exist_ok=True)

        output = mcsa_dir / "mcsa_cofactors.csv"
        return self.copy_bundled_file("mcsa_cofactors.csv", output)

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
