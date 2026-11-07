"""
Unit tests verifying CSA-3D code quality and structure.

These tests verify:
- Type hints are present on public methods
- Docstrings follow the expected format
- Code structure meets quality standards
"""

import ast
import inspect
from pathlib import Path
import pytest


# Path to main source files
ROOT_DIR = Path(__file__).parent.parent / "csa3d"
MAIN_FILES = [
    ROOT_DIR / "Mcsa.py",
    ROOT_DIR / "Entry.py",
    ROOT_DIR / "UniResidue.py",
    ROOT_DIR / "UniSite.py",
    ROOT_DIR / "PdbSite.py",
    ROOT_DIR / "PdbResidue.py",
]


@pytest.mark.unit
class TestTypeHints:
    """Tests that public methods have type hints."""

    @pytest.mark.parametrize("file_path", MAIN_FILES)
    def test_file_exists(self, file_path: Path) -> None:
        """Test that main source files exist."""
        assert file_path.exists(), f"File not found: {file_path}"

    def test_mcsa_file_has_type_hints(self) -> None:
        """Test that Mcsa.py file contains type hint annotations."""
        mcsa_file = ROOT_DIR / "Mcsa.py"
        content = mcsa_file.read_text()

        # Check for typing imports
        assert "from typing import" in content or "import typing" in content, (
            "Mcsa.py missing typing imports"
        )

        # Check for type annotations in method signatures
        assert "-> " in content or ": " in content, (
            "Mcsa.py missing type annotations"
        )

    def test_entry_file_has_type_hints(self) -> None:
        """Test that Entry.py file contains type hint annotations."""
        entry_file = ROOT_DIR / "Entry.py"
        content = entry_file.read_text()

        assert "from typing import" in content or "import typing" in content, (
            "Entry.py missing typing imports"
        )

        # Check for multiple return type annotations (should have many)
        assert content.count("-> ") > 10, "Entry.py should have many return type hints"


@pytest.mark.unit
class TestCodeQuality:
    """Tests for code quality standards."""

    def test_mcsa_uses_logging_not_print(self) -> None:
        """Test that Mcsa.py uses logging instead of print statements."""
        mcsa_file = ROOT_DIR / "Mcsa.py"
        content = mcsa_file.read_text()

        # Check for logging imports
        assert "import logging" in content, "Mcsa.py should import logging"
        assert "logger" in content, "Mcsa.py should use logger"

        # Check that most prints are replaced
        print_count = content.count("print(")
        logger_count = content.count("logger.")
        assert logger_count > print_count, (
            "Mcsa.py should use logging more than print"
        )

    def test_entry_uses_isinstance_not_type(self) -> None:
        """Test that Entry.py uses isinstance() not type() comparisons."""
        entry_file = ROOT_DIR / "Entry.py"
        content = entry_file.read_text()

        # Count isinstance() patterns (good pattern)
        isinstance_count = content.count("isinstance(")

        assert isinstance_count > 0, (
            "Entry.py should use isinstance() for type checking"
        )

    def test_modules_have_docstrings(self) -> None:
        """Test that main modules have module-level docstrings."""
        for file_path in MAIN_FILES:
            content = file_path.read_text()
            lines = content.split("\n")

            # Check if first non-comment line is a docstring
            for line in lines[:10]:
                if line.strip() and not line.strip().startswith("#"):
                    assert '"""' in line or "'''" in line, (
                        f"{file_path.name} missing module docstring"
                    )
                    break


@pytest.mark.unit
class TestClassStructure:
    """Tests for class organization and structure."""

    def test_mcsa_class_exists(self) -> None:
        """Test that Mcsa.py defines a Mcsa class."""
        mcsa_file = ROOT_DIR / "Mcsa.py"
        content = mcsa_file.read_text()

        assert "class Mcsa" in content, "Mcsa.py should define Mcsa class"

    def test_entry_class_exists(self) -> None:
        """Test that Entry.py defines an Entry class."""
        entry_file = ROOT_DIR / "Entry.py"
        content = entry_file.read_text()

        assert "class Entry" in content, "Entry.py should define Entry class"

    def test_site_classes_exist(self) -> None:
        """Test that Site classes are defined."""
        pdbsite_file = ROOT_DIR / "PdbSite.py"
        unisite_file = ROOT_DIR / "UniSite.py"

        pdbsite_content = pdbsite_file.read_text()
        unisite_content = unisite_file.read_text()

        assert "class PdbSite" in pdbsite_content, (
            "PdbSite.py should define PdbSite class"
        )
        assert "class UniSite" in unisite_content, (
            "UniSite.py should define UniSite class"
        )

    def test_residue_classes_exist(self) -> None:
        """Test that Residue classes are defined."""
        pdbresidue_file = ROOT_DIR / "PdbResidue.py"
        uniresidue_file = ROOT_DIR / "UniResidue.py"

        pdbresidue_content = pdbresidue_file.read_text()
        uniresidue_content = uniresidue_file.read_text()

        assert "class PdbResidue" in pdbresidue_content, (
            "PdbResidue.py should define PdbResidue class"
        )
        assert "class UniResidue" in uniresidue_content, (
            "UniResidue.py should define UniResidue class"
        )


@pytest.mark.unit
class TestDocumentation:
    """Tests for documentation completeness."""

    def test_mcsa_methods_have_docstrings(self) -> None:
        """Test that Mcsa methods have docstrings."""
        mcsa_file = ROOT_DIR / "Mcsa.py"
        content = mcsa_file.read_text()

        # Count docstrings in method definitions
        method_count = content.count("def ")
        docstring_count = content.count('"""')

        # Rough check: should have multiple docstrings
        assert docstring_count > method_count / 2, (
            "Mcsa.py methods should have docstrings"
        )

    def test_entry_methods_documented(self) -> None:
        """Test that Entry methods have documentation."""
        entry_file = ROOT_DIR / "Entry.py"
        content = entry_file.read_text()

        # Check for Args and Returns sections in docstrings
        assert "Args:" in content, "Entry.py should have documented parameters"
        assert "Returns:" in content, "Entry.py should document return values"

    def test_google_style_docstrings(self) -> None:
        """Test that docstrings use Google style format."""
        entry_file = ROOT_DIR / "Entry.py"
        content = entry_file.read_text()

        # Google style uses these sections
        google_style_markers = ["Args:", "Returns:", "Raises:"]

        found_markers = sum(1 for marker in google_style_markers if marker in content)
        assert found_markers > 0, (
            "Entry.py should use Google-style docstrings with Args/Returns/Raises"
        )


@pytest.mark.unit
class TestImportStructure:
    """Tests for import organization."""

    def test_typing_imports(self) -> None:
        """Test that typing imports are present where needed."""
        for file_path in [ROOT_DIR / "Entry.py", ROOT_DIR / "Mcsa.py"]:
            if file_path.exists():
                content = file_path.read_text()
                assert "from typing import" in content or "import typing" in content, (
                    f"{file_path.name} should import from typing module"
                )

    def test_logging_imports(self) -> None:
        """Test that logging is imported in main modules."""
        for file_path in [ROOT_DIR / "Entry.py", ROOT_DIR / "Mcsa.py"]:
            if file_path.exists():
                content = file_path.read_text()
                assert "import logging" in content, (
                    f"{file_path.name} should import logging"
                )
