from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional


@dataclass
class PartitionEntry:
    name: str = ""
    start: int = 1
    end: int = 0
    model: str = "AUTO"
    weight: float = 1.0


@dataclass
class PartitionFile:
    entries: List[PartitionEntry] = field(default_factory=list)

    def write_nexus(self, path: str) -> None:
        lines = ["#nexus", "begin sets;"]
        for entry in self.entries:
            lines.append(f"  charset {entry.name} = {entry.start}-{entry.end};")
        lines.append("end;")
        target = Path(path).resolve()
        target.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


@dataclass
class ConcatenatedAlignment:
    alignment: Optional[object] = None
    partition: Optional[PartitionFile] = None
    marker_order: List[str] = field(default_factory=list)
    species_order: List[str] = field(default_factory=list)
    total_sites: int = 0
    missing_data_proportion: float = 0.0


@dataclass
class AlignmentResult:
    concatenated_alignment: Optional[object] = None
    partition_file: Optional[PartitionFile] = None
    marker_order: List[str] = field(default_factory=list)
    n_markers_aligned: int = 0
    n_markers_excluded: int = 0
