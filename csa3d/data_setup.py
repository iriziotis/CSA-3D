"""Download and prepare the reference data CSA-3D needs.

Run once per data directory:

    >>> import csa3d
    >>> csa3d.setup_data()                    # defaults to ~/.csa3d/data
    >>> csa3d.setup_data("/path/to/data")     # or a location of your choice

Files that are already present are skipped, so re-running only fetches what is
missing (pass ``force=True`` to refresh everything). Assembly structures are not
downloaded here: ``Mcsa.build()`` fetches the ones it needs on demand.
"""

import csv
import gzip
import json
import logging
import shutil
import sys
import tarfile
import tempfile
import urllib.error
import urllib.request
from collections import defaultdict
from pathlib import Path
from typing import Dict, Optional, Union

from .config.settings import set_data_dir

logger = logging.getLogger(__name__)

BUNDLED_DATA_DIR = Path(__file__).parent / "bundled_data"

SIFTS_URL = "https://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/csv/{}.csv.gz"
RESOLUTION_URL = "https://files.wwpdb.org/pub/pdb/derived_data/index/resolu.idx"
RHEA_RD_URL = "https://ftp.expasy.org/databases/rhea/ctfiles/rhea-rd.tar.gz"
RHEA_TSV_URL = "https://ftp.expasy.org/databases/rhea/tsv/{}"
MCSA_API_URL = "https://www.ebi.ac.uk/thornton-srv/m-csa/api/"
PDBE_ASSEMBLY_URL = "https://www.ebi.ac.uk/pdbe/static/entry/download/{}-assembly-{}.cif.gz"


def _looks_like_html(path: Path) -> bool:
    with open(path, "rb") as f:
        head = f.read(512).lstrip().lower()
    return head.startswith(b"<!doctype html") or head.startswith(b"<html")


def _is_present(path: Path) -> bool:
    return path.is_file() and path.stat().st_size > 0 and not _looks_like_html(path)


def download_file(url: str, output_path: Path, decompress: bool = False) -> bool:
    """Download ``url`` to ``output_path``, optionally gunzipping it.

    Writes to a temporary file first and rejects HTML error pages, so a failed
    download never leaves a corrupt file behind.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    part = output_path.with_name(output_path.name + ".part")
    try:
        logger.info(f"Downloading {url}")
        with urllib.request.urlopen(url, timeout=120) as response, open(part, "wb") as out:
            source = gzip.GzipFile(fileobj=response) if decompress else response
            shutil.copyfileobj(source, out)
        if not _is_present(part):
            raise ValueError("server returned an empty file or an HTML page instead of data")
        part.replace(output_path)
        return True
    except (urllib.error.URLError, OSError, ValueError) as e:
        logger.warning(f"Failed to download {url}: {e}")
        part.unlink(missing_ok=True)
        return False


def fetch_assembly(pdb_id: str, assembly: int, dest_dir: Union[str, Path]) -> Optional[Path]:
    """Return the PDBe assembly mmCIF for ``pdb_id``, downloading it if needed."""
    output = Path(dest_dir) / f"{pdb_id}-assembly-{assembly}.cif"
    if _is_present(output):
        return output
    if download_file(PDBE_ASSEMBLY_URL.format(pdb_id, assembly), output, decompress=True):
        return output
    return None


class DataDownloader:
    """Populate a CSA-3D data directory."""

    def __init__(self, data_dir: Union[str, Path], force: bool = False):
        self.data_dir = Path(data_dir).expanduser()
        self.data_dir.mkdir(parents=True, exist_ok=True)
        self.force = force

    def _needed(self, path: Path) -> bool:
        if self.force or not _is_present(path):
            return True
        logger.info(f"Already present, skipping: {path.relative_to(self.data_dir)}")
        return False

    def _fetch(self, url: str, rel_path: str, decompress: bool = False) -> bool:
        path = self.data_dir / rel_path
        return not self._needed(path) or download_file(url, path, decompress)

    def _copy_bundled(self, filename: str, rel_path: str) -> bool:
        path = self.data_dir / rel_path
        if not self._needed(path):
            return True
        path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(BUNDLED_DATA_DIR / filename, path)
        logger.info(f"Copied bundled {filename}")
        return True

    def setup_sifts(self) -> bool:
        """SIFTS residue-level PDB <-> UniProt <-> EC mappings."""
        return all([
            self._fetch(SIFTS_URL.format(name), f"sifts/{name}.csv", decompress=True)
            for name in ("uniprot_pdb", "pdb_chain_uniprot", "pdb_chain_enzyme")
        ])

    def setup_pdbe(self) -> bool:
        """Resolution index (wwPDB) and bundled ligand/ion classification tables.

        The PDBe ion, cofactor and crystallisation-artefact tables are no longer
        published, so snapshots ship with the package.
        """
        bundled = ("pdb_ions.csv", "pdb_cofactors.csv", "crystallization_hets.csv", "metals.csv")
        return all(
            [self._fetch(RESOLUTION_URL, "pdbe/resolu.idx")]
            + [self._copy_bundled(name, f"pdbe/{name}") for name in bundled]
        )

    def setup_mcsa(self) -> bool:
        """M-CSA catalytic residues, entry metadata and (bundled) cofactors."""
        return all([
            self._fetch(MCSA_API_URL + "homologues_residues.json",
                        "mcsa/catalytic_residues_homologues.json"),
            self._fetch(MCSA_API_URL + "residues/?format=json",
                        "mcsa/catalytic_residues_references.json"),
            self._fetch_mcsa_entry_info(),
            self._copy_bundled("mcsa_cofactors.csv", "mcsa/mcsa_cofactors.csv"),
        ])

    def _fetch_mcsa_entry_info(self) -> bool:
        entry_dir = self.data_dir / "mcsa" / "entry_info"
        if not self.force and any(entry_dir.glob("mcsa_entry_info.*.json")):
            logger.info("Already present, skipping: mcsa/entry_info/")
            return True
        page = 1
        url: Optional[str] = MCSA_API_URL + "entries/?format=json&page=1"
        while url:
            path = entry_dir / f"mcsa_entry_info.{page:02d}.json"
            if not download_file(url, path):
                return False
            with open(path) as f:
                url = json.load(f).get("next")
            page += 1
        return True

    def setup_rhea(self) -> bool:
        """RHEA EC -> ChEBI reaction components, compiled from the latest release.

        Falls back to a bundled snapshot if RHEA cannot be reached.
        """
        output = self.data_dir / "rhea" / "rhea_reaction_components.csv"
        if not self._needed(output):
            return True
        with tempfile.TemporaryDirectory() as tmp:
            if _compile_rhea(Path(tmp), output):
                return True
        logger.warning("Using bundled RHEA snapshot instead")
        return self._copy_bundled("rhea_reaction_components.csv", "rhea/rhea_reaction_components.csv")

    def setup_kegg(self) -> bool:
        """KEGG EC -> compound mapping from the KEGG REST API.

        Falls back to a bundled snapshot if KEGG cannot be reached.
        """
        output = self.data_dir / "kegg" / "kegg_reaction_components.csv"
        if not self._needed(output):
            return True
        if _compile_kegg(output):
            return True
        logger.warning("Using bundled KEGG snapshot instead")
        return self._copy_bundled("kegg_reaction_components.csv", "kegg/kegg_reaction_components.csv")

    def setup_all(self) -> Dict[str, bool]:
        """Run every setup step and return {source: success}."""
        logger.info(f"Setting up CSA-3D data in {self.data_dir}")
        results = {
            "sifts": self.setup_sifts(),
            "pdbe": self.setup_pdbe(),
            "mcsa": self.setup_mcsa(),
            "rhea": self.setup_rhea(),
            "kegg": self.setup_kegg(),
        }
        for source, ok in results.items():
            logger.info(f"  {source.upper():6s} {'ok' if ok else 'FAILED'}")
        return results


def _compile_rhea(tmp: Path, output: Path) -> bool:
    rd_tar, rhea2ec, directions = tmp / "rhea-rd.tar.gz", tmp / "rhea2ec.tsv", tmp / "rhea-directions.tsv"
    if not (download_file(RHEA_RD_URL, rd_tar)
            and download_file(RHEA_TSV_URL.format("rhea2ec.tsv"), rhea2ec)
            and download_file(RHEA_TSV_URL.format("rhea-directions.tsv"), directions)):
        return False
    try:
        with tarfile.open(rd_tar, "r:gz") as tar:
            if hasattr(tarfile, "data_filter"):
                tar.extractall(tmp, filter="data")
            else:
                tar.extractall(tmp)

        directional2master = {}
        with open(directions) as f:
            next(f)
            for row in csv.reader(f, delimiter="\t"):
                directional2master[row[1]] = row[0]
                directional2master[row[2]] = row[0]

        master2ec = defaultdict(set)
        with open(rhea2ec) as f:
            next(f)
            for row in csv.reader(f, delimiter="\t"):
                master2ec[row[2]].add(row[3])

        reactions = defaultdict(set)
        for rd_file in (tmp / "rd").glob("*.rd"):
            with open(rd_file) as f:
                for line in f:
                    if line.startswith("CHEBI:"):
                        reactions[rd_file.stem].add(line.strip().split(":")[-1])

        output.parent.mkdir(parents=True, exist_ok=True)
        with open(output, "w") as f:
            f.write("RHEA_MASTER_ID,RHEA_DIRECTIONAL_ID,EC,CHEBI_IDS\n")
            # Numeric order matters: the loader keeps the last reaction per EC.
            for rhea_id in sorted(reactions, key=int):
                master = directional2master.get(rhea_id, str(int(rhea_id) - 1))
                for ec in sorted(master2ec.get(master, ())):
                    f.write(f"{master},{rhea_id},{ec},{';'.join(sorted(reactions[rhea_id]))}\n")
        logger.info(f"Compiled {len(reactions)} RHEA reactions")
        return True
    except (OSError, tarfile.TarError, ValueError, IndexError) as e:
        logger.warning(f"Could not compile RHEA data: {e}")
        return False


def _compile_kegg(output: Path) -> bool:
    from Bio.KEGG import REST

    try:
        lines = REST.kegg_link("compound", "enzyme").readlines()
    except (urllib.error.URLError, OSError) as e:
        logger.warning(f"Could not query KEGG: {e}")
        return False

    ec_compounds = defaultdict(set)
    for line in lines:
        parts = line.strip().split("\t")
        if len(parts) == 2:
            ec_compounds[parts[0].split(":")[-1]].add(parts[1].split(":")[-1])

    output.parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w") as f:
        f.write("EC,KEGG_COMPOUNDS\n")
        for ec in sorted(ec_compounds):
            f.write(f"{ec},{';'.join(sorted(ec_compounds[ec]))}\n")
    logger.info(f"Compiled {len(ec_compounds)} KEGG EC entries")
    return True


def setup_data(data_dir: Optional[Union[str, Path]] = None, force: bool = False) -> bool:
    """Download and prepare CSA-3D data, then make it the active data directory.

    Args:
        data_dir: Where to store the data. Defaults to ``~/.csa3d/data``, which
            CSA-3D finds automatically in later sessions. For any other location,
            call ``csa3d.set_data_dir(path)`` or set ``CSA3D_DATA_DIR`` in new sessions.
        force: Re-download files even if they are already present.

    Returns:
        True if every source was set up successfully.
    """
    data_dir = Path(data_dir).expanduser() if data_dir else Path.home() / ".csa3d" / "data"
    results = DataDownloader(data_dir, force=force).setup_all()
    set_data_dir(str(data_dir))
    failed = [source for source, ok in results.items() if not ok]
    if failed:
        logger.warning(f"Setup incomplete; failed sources: {', '.join(failed)}. "
                       "Re-run setup_data() to retry only the missing files.")
    return not failed


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
    sys.exit(0 if setup_data(sys.argv[1] if len(sys.argv) > 1 else None) else 1)
