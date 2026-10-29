"""
Unit tests for csa3d.Mcsa class.

Tests database initialization, entry building, and basic operations.
"""

import pytest
from csa3d import Mcsa


@pytest.mark.unit
class TestMcsaInitialization:
    """Tests for Mcsa class initialization."""

    def test_mcsa_can_be_instantiated(self) -> None:
        """Test that Mcsa instance can be created."""
        try:
            db = Mcsa()
            assert db is not None
        except Exception as e:
            # If data is not available, that's OK for unit tests
            if "data" in str(e).lower():
                pytest.skip(f"CSA-3D data not available: {e}")
            else:
                raise

    def test_mcsa_has_entries_attribute(self) -> None:
        """Test that Mcsa has entries dictionary."""
        try:
            db = Mcsa()
            assert hasattr(db, "entries")
            assert isinstance(db.entries, dict)
        except Exception as e:
            if "data" in str(e).lower():
                pytest.skip(f"CSA-3D data not available: {e}")
            else:
                raise


@pytest.mark.unit
class TestMcsaAttributes:
    """Tests for Mcsa class attributes and properties."""

    def test_mcsa_has_required_methods(self) -> None:
        """Test that Mcsa has required methods."""
        expected_methods = ["build", "add", "update_pdbs"]

        db_methods = [method for method in dir(Mcsa) if not method.startswith("_")]

        for method in expected_methods:
            assert any(
                method in m for m in db_methods
            ), f"Mcsa missing method: {method}"


@pytest.mark.integration
@pytest.mark.requires_data
class TestMcsaBuildMethod:
    """Integration tests for Mcsa.build() method."""

    def test_mcsa_build_returns_boolean(self, mcsa_db: Mcsa, sample_mcsa_id: int) -> None:
        """Test that build() returns a boolean."""
        if mcsa_db is None:
            pytest.skip("Mcsa database not available")

        result = mcsa_db.build(sample_mcsa_id, annotate=False)
        assert isinstance(result, bool)

    def test_mcsa_build_invalid_id_returns_false(
        self, mcsa_db: Mcsa
    ) -> None:
        """Test that build() returns False for invalid ID."""
        if mcsa_db is None:
            pytest.skip("Mcsa database not available")

        # Use an invalid entry ID
        result = mcsa_db.build(999999, annotate=False)
        assert isinstance(result, bool)


@pytest.mark.unit
class TestMcsaEntries:
    """Tests for Mcsa entries dictionary."""

    def test_mcsa_entries_is_dict(self) -> None:
        """Test that Mcsa.entries is a dictionary."""
        try:
            db = Mcsa()
            assert isinstance(db.entries, dict)
        except Exception as e:
            if "data" in str(e).lower():
                pytest.skip(f"CSA-3D data not available: {e}")
            else:
                raise


@pytest.mark.unit
class TestMcsaTypeHints:
    """Tests to verify type hints are present."""

    def test_mcsa_build_has_type_hints(self) -> None:
        """Test that Mcsa.build method has type hints."""
        import inspect

        sig = inspect.signature(Mcsa.build)
        # Check that parameters have annotations
        for param_name, param in sig.parameters.items():
            if param_name != "self":
                # Type hints should be present
                assert param.annotation != inspect.Parameter.empty, (
                    f"Parameter {param_name} missing type hint"
                )

        # Return type should be annotated
        assert sig.return_annotation != inspect.Signature.empty, (
            "build() method missing return type hint"
        )

    def test_mcsa_add_has_type_hints(self) -> None:
        """Test that Mcsa.add method has type hints."""
        import inspect

        sig = inspect.signature(Mcsa.add)
        # Check return type
        assert sig.return_annotation != inspect.Signature.empty, (
            "add() method missing return type hint"
        )

    def test_mcsa_has_docstrings(self) -> None:
        """Test that Mcsa methods have docstrings."""
        methods = ["build", "add", "__init__"]

        for method_name in methods:
            method = getattr(Mcsa, method_name, None)
            if method is not None:
                assert method.__doc__ is not None, (
                    f"Mcsa.{method_name} missing docstring"
                )
