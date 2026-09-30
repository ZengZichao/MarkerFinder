"""Unit tests for database version management."""

from pathlib import Path

from markerfinder.utils.db_versioning import DatabaseVersionManager


class TestDatabaseVersionManager:
    def test_compute_file_hash(self, tmp_path):
        path = tmp_path / "data.txt"
        path.write_text("hello")
        mgr = DatabaseVersionManager(str(tmp_path))
        h1 = mgr.compute_file_hash(str(path))
        h2 = mgr.compute_file_hash(str(path))
        assert len(h1) == 64
        assert h1 == h2

    def test_verify_missing_database(self, tmp_path):
        mgr = DatabaseVersionManager(str(tmp_path))
        info = mgr.verify_database("missing", str(tmp_path / "no_such_dir"))
        assert info["exists"] is False
        assert info["hash"] == ""

    def test_verify_directory_hash(self, tmp_path):
        db_dir = tmp_path / "gtdb_markers"
        db_dir.mkdir()
        (db_dir / "a.faa").write_text("seq1")
        (db_dir / "b.faa").write_text("seq2")
        mgr = DatabaseVersionManager(str(tmp_path))
        info = mgr.verify_database("gtdb_markers", str(db_dir))
        assert info["exists"] is True
        assert len(info["hash"]) == 64

    def test_verify_file_hash(self, tmp_path):
        db_file = tmp_path / "db.bin"
        db_file.write_bytes(b"data")
        mgr = DatabaseVersionManager(str(tmp_path))
        info = mgr.verify_database("db", str(db_file))
        assert info["exists"] is True
        assert len(info["hash"]) == 64

    def test_generate_version_report(self, tmp_path):
        mgr = DatabaseVersionManager(str(tmp_path))
        versions = {"gtdb": {"name": "gtdb", "hash": "abc", "exists": True}}
        report = mgr.generate_version_report(versions)
        assert report["database_versions"] == versions
        assert "markerfinder_version" in report

    def test_get_gtdb_mapping_when_readme_exists(self, tmp_path):
        gtdb_dir = tmp_path / "db" / "gtdb_markers"
        gtdb_dir.mkdir(parents=True)
        (gtdb_dir / "README.md").write_text("mapping")
        mgr = DatabaseVersionManager(str(tmp_path / "db"))
        mapping = mgr._get_gtdb_mapping()
        assert mapping["exists"] is True

    def test_get_gtdb_mapping_when_readme_missing(self, tmp_path):
        mgr = DatabaseVersionManager(str(tmp_path / "db"))
        mapping = mgr._get_gtdb_mapping()
        assert mapping["exists"] is False
