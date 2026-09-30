import pytest

from markerfinder.models.alignment import (
    ConcatenatedAlignment,
    PartitionEntry,
    PartitionFile,
    AlignmentResult,
)


class TestPartitionEntry:
    def test_creation(self):
        entry = PartitionEntry(name="COG001", start=1, end=300, model="LG+F")
        assert entry.name == "COG001"
        assert entry.start == 1
        assert entry.end == 300
        assert entry.model == "LG+F"

    def test_default_model(self):
        entry = PartitionEntry(name="test", start=1, end=100)
        assert entry.model == "AUTO"


class TestPartitionFile:
    def test_write_nexus(self, tmp_path):
        entries = [
            PartitionEntry(name="COG001", start=1, end=300),
            PartitionEntry(name="COG002", start=301, end=600),
        ]
        pf = PartitionFile(entries=entries)
        out = tmp_path / "test.nex"
        pf.write_nexus(str(out))
        content = out.read_text()
        assert "#nexus" in content
        assert "charset COG001 = 1-300;" in content
        assert "charset COG002 = 301-600;" in content


class TestConcatenatedAlignment:
    def test_default(self):
        ca = ConcatenatedAlignment()
        assert ca.total_sites == 0
        assert ca.missing_data_proportion == 0.0
        assert ca.marker_order == []


class TestAlignmentResult:
    def test_default(self):
        ar = AlignmentResult()
        assert ar.n_markers_aligned == 0
        assert ar.n_markers_excluded == 0
