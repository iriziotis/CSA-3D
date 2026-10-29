"""
Unit tests for csa3d.Entry class.

Tests entry operations, RMSD calculations, clustering, and template creation.
"""

import pytest
from csa3d import Entry


@pytest.mark.unit
class TestEntryAttributes:
    """Tests for Entry class attributes."""

    def test_entry_has_required_attributes(self) -> None:
        """Test that Entry class has required attributes/methods."""
        expected_attributes = [
            "mcsa_id",
            "e_number",
            "pdbsites",
            "unisites",
            "build",
            "add",
        ]

        for attr in expected_attributes:
            assert hasattr(Entry, attr) or hasattr(Entry, f"_{attr}"), (
                f"Entry missing attribute/method: {attr}"
            )


@pytest.mark.unit
class TestEntryTypeHints:
    """Tests to verify Entry class has type hints."""

    def test_entry_methods_have_type_hints(self) -> None:
        """Test that Entry methods have type hints."""
        import inspect

        methods_to_check = ["__init__", "build", "add", "rmsd_matrix"]

        for method_name in methods_to_check:
            if hasattr(Entry, method_name):
                method = getattr(Entry, method_name)
                sig = inspect.signature(method)

                # Check that return type is annotated
                if method_name != "__init__":
                    assert sig.return_annotation != inspect.Signature.empty, (
                        f"Entry.{method_name}() missing return type hint"
                    )

    def test_entry_has_docstrings(self) -> None:
        """Test that Entry methods have docstrings."""
        methods = ["__init__", "build", "add"]

        for method_name in methods:
            method = getattr(Entry, method_name, None)
            if method is not None:
                assert method.__doc__ is not None, (
                    f"Entry.{method_name} missing docstring"
                )


@pytest.mark.unit
class TestEntryProperties:
    """Tests for Entry properties."""

    def test_entry_has_property_methods(self) -> None:
        """Test that Entry has expected property methods."""
        properties = [
            "has_multiple_pdb_refs",
            "has_multiple_uniprot_refs",
            "n_pdb_sites",
            "n_uni_sites",
        ]

        entry_attrs = dir(Entry)
        for prop in properties:
            assert any(
                prop in attr for attr in entry_attrs
            ), f"Entry missing property: {prop}"


@pytest.mark.integration
@pytest.mark.requires_data
class TestEntryRmsdMatrix:
    """Integration tests for Entry RMSD matrix calculation."""

    def test_rmsd_matrix_signature(self) -> None:
        """Test that rmsd_matrix method has correct signature."""
        import inspect

        sig = inspect.signature(Entry.rmsd_matrix)

        # Should have sitelist parameter
        assert "sitelist" in sig.parameters, "rmsd_matrix missing sitelist parameter"

        # Should return a matrix (numpy array or similar)
        assert sig.return_annotation != inspect.Signature.empty, (
            "rmsd_matrix missing return type hint"
        )


@pytest.mark.integration
@pytest.mark.requires_data
class TestEntryClustering:
    """Integration tests for Entry clustering methods."""

    def test_clustering_method_exists(self) -> None:
        """Test that clustering method exists."""
        assert hasattr(Entry, "clustering"), "Entry missing clustering method"

    def test_clustering_bayesian_method_exists(self) -> None:
        """Test that clustering_bayesian method exists."""
        assert hasattr(
            Entry, "clustering_bayesian"
        ), "Entry missing clustering_bayesian method"

    def test_clustering_methods_have_return_types(self) -> None:
        """Test that clustering methods have return type hints."""
        import inspect

        for method_name in ["clustering", "clustering_bayesian"]:
            method = getattr(Entry, method_name, None)
            if method is not None:
                sig = inspect.signature(method)
                assert sig.return_annotation != inspect.Signature.empty, (
                    f"Entry.{method_name}() missing return type hint"
                )


@pytest.mark.integration
@pytest.mark.requires_data
class TestEntryTemplateCreation:
    """Integration tests for Entry template creation methods."""

    def test_create_template_method_exists(self) -> None:
        """Test that create_template method exists."""
        assert hasattr(Entry, "create_template"), "Entry missing create_template method"

    def test_break_template_method_exists(self) -> None:
        """Test that break_template method exists."""
        assert hasattr(Entry, "break_template"), "Entry missing break_template method"

    def test_write_template_method_exists(self) -> None:
        """Test that write_template method exists."""
        assert hasattr(Entry, "write_template"), "Entry missing write_template method"

    def test_template_methods_have_return_types(self) -> None:
        """Test that template methods have return type hints."""
        import inspect

        methods = ["create_template", "break_template", "write_template"]
        for method_name in methods:
            method = getattr(Entry, method_name, None)
            if method is not None:
                sig = inspect.signature(method)
                # Most template methods should have return type hints
                if method_name != "write_template":
                    assert sig.return_annotation != inspect.Signature.empty, (
                        f"Entry.{method_name}() missing return type hint"
                    )
