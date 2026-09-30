"""Database version management utilities."""

import hashlib
import json
import logging
from pathlib import Path
from typing import Dict, Optional

from markerfinder._version import __version__


logger = logging.getLogger(__name__)


class DatabaseVersionManager:
    """数据库版本管理器

    通过 SHA256 哈希校验锁定数据库版本，确保分析可复现。
    """

    def __init__(self, db_dir: str):
        self.db_dir = Path(db_dir)

    def compute_file_hash(self, file_path: str) -> str:
        """计算文件的 SHA256 哈希"""
        sha256 = hashlib.sha256()
        with open(file_path, "rb") as f:
            for chunk in iter(lambda: f.read(8192), b""):
                sha256.update(chunk)
        return sha256.hexdigest()

    # Optional per-database expected hashes, loaded from
    # Db/expected_hashes.json ({"<db_name>": "<sha256>",...}).
    EXPECTED_HASHES_FILE = "expected_hashes.json"

    def _load_expected_hashes(self) -> Dict:
        import json as _json

        expected_file = self.db_dir / self.EXPECTED_HASHES_FILE
        if not expected_file.exists():
            return {}
        try:
            return _json.loads(expected_file.read_text(encoding="utf-8")) or {}
        except Exception:
            return {}

    def verify_database(self, db_name: str, db_path: str) -> Dict:
        """验证数据库并返回版本信息。

        when ``db/expected_hashes.json`` carries an expected
        hash for ``db_name``, the computed hash is COMPARED against it and the
        result carries ``expected_hash`` / ``hash_match`` — a database that
        "computes" but no longer "matches" is reported, never silently passed.

        Returns:
            {'name', 'path', 'hash', 'exists',
             'expected_hash': str|None, 'hash_match': bool|None}
        """
        path = Path(db_path)
        if not path.exists():
            return {
                "name": db_name,
                "path": db_path,
                "hash": "",
                "exists": False,
                "expected_hash": self._load_expected_hashes().get(db_name),
                "hash_match": None,
            }

        if path.is_dir():
            # 对目录，计算所有文件的组合哈希
            all_hashes = []
            for f in sorted(path.rglob("*")):
                if f.is_file():
                    all_hashes.append(self.compute_file_hash(str(f)))
            combined = hashlib.sha256("".join(all_hashes).encode()).hexdigest()
        else:
            combined = self.compute_file_hash(db_path)

        expected = self._load_expected_hashes().get(db_name)
        hash_match = (combined == expected) if expected else None
        return {
            "name": db_name,
            "path": db_path,
            "hash": combined,
            "exists": True,
            "expected_hash": expected,
            "hash_match": hash_match,
        }

    def generate_version_report(self, db_versions: Dict[str, Dict]) -> Dict:
        """生成数据库版本报告"""
        return {
            "database_versions": db_versions,
            "gtdb_mapping": self._get_gtdb_mapping(),
            "markerfinder_version": __version__,
        }

    def _get_gtdb_mapping(self) -> Dict:
        """获取 GTDB 标记集映射信息"""
        mapping_file = self.db_dir / "gtdb_markers" / "README.md"
        if mapping_file.exists():
            return {"mapping_file": str(mapping_file), "exists": True}
        return {"mapping_file": str(mapping_file), "exists": False}
