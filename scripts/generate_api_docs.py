#!/usr/bin/env python3
"""
Generate API reference documentation from docstrings.

Scans Python source files and extracts class and method documentation
to generate an API reference markdown file.
"""

import re
import ast
from pathlib import Path
from typing import List, Dict, Tuple, Optional


class DocExtractor:
    """Extract docstrings from Python files."""

    def __init__(self, file_path: Path):
        """Initialize extractor with a Python file."""
        self.file_path = file_path
        self.content = file_path.read_text()
        self.tree = ast.parse(self.content)

    def extract_module_docstring(self) -> Optional[str]:
        """Extract module-level docstring."""
        return ast.get_docstring(self.tree)

    def extract_classes(self) -> List[Tuple[str, str, List[Tuple[str, str]]]]:
        """Extract classes and their methods with docstrings.

        Returns:
            List of (class_name, class_doc, methods) tuples
            where methods is list of (method_name, method_doc)
        """
        classes = []

        for node in ast.walk(self.tree):
            if isinstance(node, ast.ClassDef):
                class_doc = ast.get_docstring(node) or "No documentation"
                methods = []

                for item in node.body:
                    if isinstance(item, ast.FunctionDef) and not item.name.startswith("_"):
                        method_doc = ast.get_docstring(item) or "No documentation"
                        methods.append((item.name, method_doc))

                if methods:
                    classes.append((node.name, class_doc, methods))

        return classes


def generate_api_docs(source_dir: Path, output_file: Path) -> None:
    """Generate API documentation from Python files.

    Args:
        source_dir: Directory containing Python source files
        output_file: Path to output markdown file
    """
    # Find main source files
    main_files = [
        "Mcsa.py",
        "Entry.py",
        "UniSite.py",
        "UniResidue.py",
        "PdbSite.py",
        "PdbResidue.py",
    ]

    doc_lines = [
        "# CSA-3D API Reference",
        "",
        "Auto-generated documentation from source code docstrings.",
        "",
        "---",
        "",
    ]

    for filename in main_files:
        file_path = source_dir / filename

        if not file_path.exists():
            continue

        print(f"Processing {filename}...")

        extractor = DocExtractor(file_path)

        # Add file section
        doc_lines.append(f"## {filename}")
        doc_lines.append("")

        module_doc = extractor.extract_module_docstring()
        if module_doc:
            doc_lines.append(module_doc)
            doc_lines.append("")

        # Extract and document classes
        classes = extractor.extract_classes()

        for class_name, class_doc, methods in classes:
            doc_lines.append(f"### Class: `{class_name}`")
            doc_lines.append("")
            doc_lines.append(class_doc)
            doc_lines.append("")

            if methods:
                doc_lines.append(f"**Methods ({len(methods)}):**")
                doc_lines.append("")

                for method_name, method_doc in methods:
                    doc_lines.append(f"#### `{method_name}()`")
                    doc_lines.append("")
                    doc_lines.append(method_doc)
                    doc_lines.append("")

        doc_lines.append("---")
        doc_lines.append("")

    # Write output
    output_file.write_text("\n".join(doc_lines))
    print(f"\n✓ API documentation written to {output_file}")


def main():
    """Main entry point."""
    root_dir = Path(__file__).parent
    source_dir = root_dir

    output_file = root_dir / "docs" / "api_reference.md"
    output_file.parent.mkdir(exist_ok=True)

    print("=" * 60)
    print("CSA-3D API Reference Generator")
    print("=" * 60)
    print()

    generate_api_docs(source_dir, output_file)

    print()
    print("=" * 60)
    print("✓ Documentation generation complete!")
    print("=" * 60)


if __name__ == "__main__":
    main()
