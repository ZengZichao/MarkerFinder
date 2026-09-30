"""Unit tests for custom exception hierarchy."""

import pytest

from markerfinder.exceptions import (
    PhyloToolError,
    PhyloFormatError,
    TaxonomyConflictError,
    MonophylyError,
    TreeValidationError,
    SequenceValidationError,
    CrossValidationError,
    InputError,
    ConfigError,
    MultiTreeError,
)


class TestExceptionHierarchy:
    def test_all_inherit_from_base(self):
        exceptions = [
            PhyloFormatError,
            TaxonomyConflictError,
            MonophylyError,
            TreeValidationError,
            SequenceValidationError,
            CrossValidationError,
            InputError,
            ConfigError,
            MultiTreeError,
        ]
        for exc_class in exceptions:
            assert issubclass(exc_class, PhyloToolError)

    def test_catch_with_base(self):
        with pytest.raises(PhyloToolError):
            raise TreeValidationError("test tree error")

    def test_specific_catch(self):
        with pytest.raises(MonophylyError):
            raise MonophylyError("taxon not found")

    def test_error_message(self):
        try:
            raise ConfigError("bad config")
        except ConfigError as e:
            assert "bad config" in str(e)
