"""Shared fixtures for the MarkerFinder validation (full-feature E2E) suite.

Design rules this file enforces
-------------------------------
* Every case gets its own output directory, named after the pytest test id.
  Nothing is shared between cases, so the suite can run under
  ``pytest -n N`` (and re-run single cases) without cross-talk.
* The suite never writes into ``data/``. Runtime-materialised inputs (a
  directory of symlinks to the shipped genome pool, empty scratch
  directories) live under ``.work/``, which is disposable and git-ignored.
* A case fails with a message that says what to re-run, never with a bare
  ``assert False``.
"""

from __future__ import annotations

import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence

import pytest

ROOT = Path(__file__).parent.resolve()          # Validation/
REPO = ROOT.parent
DATA = ROOT / "data"
WORK = ROOT / ".work"
OUTPUTS = WORK / "outputs"
RESULTS = ROOT / "results"

MARKERFINDER_PKG = REPO / "markerfinder"


def pytest_configure(config):
    config.addinivalue_line(
        "markers", "slow: needs the 12-genome benchmark set (deselect with -m 'not slow')"
    )
    config.addinivalue_line(
        "markers", "capability(name): the capability a case is registered against "
                   "in the coverage matrix (see validation/capabilities.py)"
    )


def pytest_collection_modifyitems(session, config, items):
    """Dump ``capability id -> test node ids`` for the coverage ratchet.

    Written per xdist worker and merged on read (see
    ``capabilities.covered_capabilities``), so the matrix reflects the tests
    that actually exist in this checkout rather than a list someone maintains
    beside them. Collecting only a subset of the suite produces a subset
    snapshot; ``test_v90`` refuses to run on a partial collection so a
    narrowed invocation can never make coverage look complete.
    """
    worker = os.environ.get("PYTEST_XDIST_WORKER", "main")
    mapping: Dict[str, List[str]] = {}
    for item in items:
        for marker in item.iter_markers(name="capability"):
            for name in marker.args:
                mapping.setdefault(str(name), [])
                if item.nodeid not in mapping[str(name)]:
                    mapping[str(name)].append(item.nodeid)
    SNAPSHOT_DIR = WORK / "capability_snapshot"
    SNAPSHOT_DIR.mkdir(parents=True, exist_ok=True)
    path = SNAPSHOT_DIR / f"{worker}.json"
    existing = {}
    if path.exists():
        try:
            existing = json.loads(path.read_text(encoding="utf-8"))
        except json.JSONDecodeError:
            existing = {}
    existing.update({k: sorted(v) for k, v in mapping.items()})
    payload = dict(sorted(existing.items()))
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    # Also expose the collection size for the "was the suite collected whole"
    # Guard in test_v90
    (SNAPSHOT_DIR / f"{worker}.count").write_text(str(len(items)), encoding="utf-8")


def _worker_dir(base: Path) -> Path:
    """Per-xdist-worker directory so parallel prep cannot race."""
    worker = os.environ.get("PYTEST_XDIST_WORKER", "main")
    path = base / worker
    path.mkdir(parents=True, exist_ok=True)
    return path


@pytest.fixture(scope="session")
def root_dir() -> Path:
    return ROOT


@pytest.fixture(scope="session")
def repo_root() -> Path:
    return REPO


@pytest.fixture(scope="session")
def data_dir() -> Path:
    if not DATA.exists():
        raise FileNotFoundError(
            f"{DATA} is missing. Build the test data first:\n"
            "  cd validation/scripts && python 01_fetch_genomes.py && "
            "python 02_build_marker_sets.py && python 03_build_fixtures.py"
        )
    return DATA


@pytest.fixture(scope="session")
def results_dir() -> Path:
    RESULTS.mkdir(parents=True, exist_ok=True)
    return RESULTS


# --------------------------------------------------------------------------
# Data locations
# --------------------------------------------------------------------------
def _require(path: Path, hint: str) -> Path:
    if not path.exists():
        raise FileNotFoundError(f"{path} missing — {hint}")
    return path


@pytest.fixture(scope="session")
def genome_pool(data_dir: Path) -> Path:
    return _require(data_dir / "genomes", "run 01_fetch_genomes.py")


@pytest.fixture(scope="session")
def provenance(data_dir: Path) -> List[Dict[str, str]]:
    path = _require(data_dir / "PROVENANCE.tsv", "run 01_fetch_genomes.py")
    lines = path.read_text(encoding="utf-8").rstrip("\n").split("\n")
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


@pytest.fixture(scope="session")
def sets_dir(data_dir: Path) -> Path:
    return _require(data_dir / "sets", "run 01_fetch_genomes.py")


@pytest.fixture(scope="session")
def taxonomy_table(data_dir: Path):
    def _get(name: str = "small8") -> Path:
        return _require(data_dir / "taxonomy" / f"taxonomy_{name}.tsv",
                        "run 03_build_fixtures.py")
    return _get


@pytest.fixture(scope="session")
def markers_dir(data_dir: Path):
    def _get(name: str = "small8_core") -> Path:
        return _require(data_dir / "markers" / name,
                        "run 02_build_marker_sets.py")
    return _get


@pytest.fixture(scope="session")
def hmm_dir(data_dir: Path):
    def _get(name: str = "core") -> Path:
        return _require(data_dir / "hmms" / name, "run 02_build_marker_sets.py")
    return _get


@pytest.fixture(scope="session")
def marker_set(data_dir: Path):
    """``marker_set("small8")`` -> the complete(core) markers of that set.

    ``flavour="spanning"`` gives the occupancy-spanning directory, whose
    markers are not all present in every genome — the only way an occupancy
    threshold can be shown to reject as well as to accept.
    """
    def _get(set_name: str, flavour: str = "core") -> Path:
        name = set_name if flavour == "spanning" else f"{set_name}_{flavour}"
        return _require(data_dir / "markers" / name,
                        "run 02_build_marker_sets.py")
    return _get


# --------------------------------------------------------------------------
# Input materialisation
# --------------------------------------------------------------------------
@pytest.fixture(scope="session")
def genome_input(sets_dir: Path, genome_pool: Path):
    """Materialise an input directory of.faa symlinks for a named set.

    ``naming='embedded'`` gives the Format A file names
    (``<acc>_d_Bacteria_p_..._g_<genus>.faa``) the embedded-taxonomy
    auto-detection reads, so that path is exercised on the same real
    sequences rather than on a second copy of the data.
    """
    cache: Dict[tuple, Path] = {}

    def _materialise(set_name: str, naming: str = "plain") -> Path:
        key = (set_name, naming)
        if key in cache:
            return cache[key]
        members = [
            x.strip() for x in (sets_dir / f"{set_name}.txt")
            .read_text(encoding="utf-8").splitlines() if x.strip()
        ]
        provenance = {}
        prot = sets_dir.parent / "PROVENANCE.tsv"
        if prot.exists() and naming == "embedded":
            lines = prot.read_text(encoding="utf-8").rstrip("\n").split("\n")
            header = lines[0].split("\t")
            for line in lines[1:]:
                row = dict(zip(header, line.split("\t")))
                provenance[row["accession"]] = row

        worker = _worker_dir(WORK / "inputs")
        target = worker / f"{set_name}_{naming}"
        if target.exists():
            shutil.rmtree(target)
        target.mkdir(parents=True)
        for acc in members:
            src = genome_pool / f"{acc}.faa"
            if not src.exists():
                raise FileNotFoundError(f"{src} missing — run 01_fetch_genomes.py")
            if naming == "plain":
                dst = target / f"{acc}.faa"
            else:
                row = provenance[acc]
                # Canonical Format A label: _d_X_p_Y_c_Z_o_W_f_V_g_U
                pairs = [("d", row["domain"]), ("p", row["phylum"]),
                         ("c", row["class"]), ("o", row["order"]),
                         ("f", row["family"]), ("g", row["genus"])]
                seg = "".join(f"_{code}_{value.replace(' ', '_')}"
                              for code, value in pairs)
                dst = target / f"{acc}{seg}.faa"
            dst.symlink_to(src)
        cache[key] = target
        return target

    return _materialise


@pytest.fixture(scope="session")
def empty_input_dir():
    """A directory that exists and contains nothing (an edge case, not a bug)."""
    def _get(name: str = "empty") -> Path:
        path = _worker_dir(WORK / "empty") / name
        path.mkdir(parents=True, exist_ok=True)
        for child in path.iterdir():
            if child.is_file():
                child.unlink()
        return path
    return _get


# --------------------------------------------------------------------------
# The CLI runner
# --------------------------------------------------------------------------
@dataclass
class Run:
    rc: int
    stdout: str
    stderr: str
    argv: List[str]
    out_dir: Optional[Path]

    @property
    def text(self) -> str:
        return self.stdout + self.stderr

    def tail(self, n: int = 60) -> str:
        lines = self.text.strip().splitlines()
        return "\n".join(lines[-n:])

    def assert_ok(self, why: str = "") -> "Run":
        assert self.rc == 0, (
            f"markerfinder exited {self.rc}. {why}\n$ markerfinder "
            f"{shlex.join(self.argv)}\n--- last output ---\n{self.tail()}"
        )
        return self

    def product(self, relative: str) -> Path:
        assert self.out_dir is not None, "run had no -o output directory"
        path = self.out_dir / relative
        if not path.exists():
            listing = sorted(
                str(p.relative_to(self.out_dir)) for p in self.out_dir.rglob("*")
                if p.is_file()
            )
            raise AssertionError(
                f"expected product {relative!r} is missing.\nFiles present: "
                f"{listing}\n--- last output ---\n{self.tail()}"
            )
        return path

    def recorded_parameters(self) -> dict:
        """The ``parameters`` object the run wrote into run_config.json.

        Keys are the pipeline's own config groups (``selection_config``,
        ``hgt_config``, ``phylo_config``, ``taxonomy_config`` …) rather than
        CLI option names, so a case must say WHERE it expects a value. Reading
        through this one accessor keeps that explicit instead of leaving each
        case to guess the nesting.
        """
        payload = json.loads(
            self.product("Phase5_metadata/run_config.json")
            .read_text(encoding="utf-8")
        )
        return payload.get("parameters", payload)

    def recorded(self, group: str, key: str):
        """``run.recorded('selection_config', 'max_markers')``.

        A top-level key is found too, so ``recorded('', 'cpus')`` works for the
        groups that are not nested.
        """
        params = self.recorded_parameters()
        if group:
            block = params.get(group)
            assert isinstance(block, dict), (
                f"run_config.json has no {group!r} block: {sorted(params)}"
            )
            return block.get(key)
        return params.get(key)


@pytest.fixture(scope="session", autouse=True)
def _environment_on_path():
    """Put the invoking interpreter's ``bin`` on PATH for the whole session."""
    bindir = str(Path(sys.executable).parent)
    previous = os.environ.get("PATH", "")
    if bindir not in previous.split(os.pathsep):
        os.environ["PATH"] = f"{bindir}{os.pathsep}{previous}" if previous else bindir
    yield
    os.environ["PATH"] = previous


def _tool_env(run_tmp: Optional[Path] = None) -> Dict[str, str]:
    """Environment for the subprocess runs: PATH, tool homes, and a writable TMPDIR.

    Three things a bare ``python -m pytest`` inherits wrongly:

    * PATH — bioconda installs mafft/trimal/IQ-TREE/ASTRAL/FastTree/hmmsearch
      into the environment's ``bin``, so a run must see it. Without this the
      same command passes in an activated shell and fails with "mafft not found
      on PATH" otherwise, which makes the suite look broken when the shell is
      merely un-activated.
    * TMPDIR — MAFFT creates its scratch directory with ``mktemp -dt`` and,
      when that fails, does NOT stop: it continues with an empty path, produces
      no alignment, and the pipeline then reports "No markers could be
      aligned". A read-only or full ``/tmp`` therefore turns into a wrong
      scientific verdict rather than a clear error. Pinning TMPDIR inside
      ``validation/.work`` makes each run's scratch space explicit, writable,
      and removable with the rest of the work directory.
    * ``MAFFT_TMPDIR`` is set to the same place so MAFFT does not depend on the
      ambient TMPDIR at all.
    """
    env = dict(os.environ)
    prefix = Path(sys.executable).parent.parent
    bindir = str(prefix / "bin")
    existing = env.get("PATH", "")
    if bindir not in existing.split(os.pathsep):
        env["PATH"] = bindir + os.pathsep + existing if existing else bindir
    env["PYTHONUNBUFFERED"] = "1"

    if run_tmp is not None:
        run_tmp.mkdir(parents=True, exist_ok=True)
        env["TMPDIR"] = str(run_tmp)
        env["MAFFT_TMPDIR"] = str(run_tmp)

    if not env.get("MAFFT_HOME"):
        for candidate in (prefix / "opt" / "mafft", prefix / "libexec" / "mafft"):
            if candidate.is_dir():
                env["MAFFT_HOME"] = str(candidate)
                break
    if not env.get("CHECKM_DATA") and (Path.home() / ".checkm").is_dir():
        env["CHECKM_DATA"] = str(Path.home() / ".checkm")
    return env


@pytest.fixture(scope="session")
def run_markerfinder(repo_root: Path):
    """Invoke the CLI as a subprocess, exactly as a user would."""
    counter = {"n": 0}

    def _run(args: Sequence[str], out_name: Optional[str] = None,
             cwd: Optional[Path] = None, expect_rc: Optional[int] = None,
             timeout: int = 3600) -> Run:
        argv = [sys.executable, "-m", "markerfinder"] + [str(a) for a in args]
        out_dir = None
        for i, token in enumerate(argv):
            if token == "-o" and i + 1 < len(argv):
                out_dir = Path(argv[i + 1])
        run_tmp = _worker_dir(WORK / "tmp")
        proc = subprocess.run(
            argv, cwd=str(cwd or repo_root), capture_output=True, text=True,
            timeout=timeout, env=_tool_env(run_tmp),
        )
        if expect_rc is not None:
            assert proc.returncode == expect_rc, (
                f"expected exit code {expect_rc}, got {proc.returncode}\n"
                f"$ {' '.join(argv)}\n--- last output ---\n"
                + "\n".join((proc.stdout + proc.stderr).strip().splitlines()[-40:])
            )
        counter["n"] += 1
        return Run(rc=proc.returncode, stdout=proc.stdout, stderr=proc.stderr,
                   argv=list(argv[2:]), out_dir=out_dir)

    return _run


@pytest.fixture(scope="session")
def mf(run_markerfinder, repo_root: Path, data_dir: Path, taxonomy_table,
       marker_set, genome_input):
    """The standard cheap run, with per-case output directories.

    Give it a set name (``mf("small8")``) and the input directory, marker
    directory and taxonomy table are taken from that set; give it a path and
    only ``-i`` is yours. Everything else has one documented default so a
    case's command line shows only what the case is testing:
    4 markers from the complete core set, FastTree gene trees, coalescent off,
    CheckM skipped — measured wall clock 1m34s for a 4-genome / 6-marker run
    against 7m19s for the same dataset over 30 markers.
    """
    counter = {"n": 0}

    def _mf(target="quad4", *, markers: Optional[Path] = None,
            taxonomy: Optional[Path] = None,
            extra: Sequence[str] = (), omit: Sequence[str] = (),
            out_dir: Optional[Path] = None, expect_rc: Optional[int] = None,
            timeout: int = 3600, **kwargs) -> Run:
        counter["n"] += 1
        if isinstance(target, Path):
            input_dir = target
            set_name = "quad4"
        else:
            input_dir = genome_input(target)
            set_name = target
        worker = os.environ.get("PYTEST_XDIST_WORKER", "main")
        out = out_dir or (OUTPUTS / worker / f"auto_{counter['n']:04d}")
        out.mkdir(parents=True, exist_ok=True)
        chosen_markers = markers or marker_set(set_name)
        chosen_taxonomy = taxonomy or taxonomy_table(set_name)
        args: List[str] = ["-i", str(input_dir), "-o", str(out), "-t", "2"]
        defaults = {
            "--marker-mode": "gtdb_tk",
            "--gtdb-markers-dir": str(chosen_markers),
            "--taxonomy-table": str(chosen_taxonomy),
            "--max-markers": "4",
            # No --monophyly-rank here on purpose: the suite's baseline run uses
            # The documented default (auto = one rank below the taxonomic scope
            # Of the tree's tips). Pinning a rank here would make the default run
            # Of a cohesive 4-genome set unmeasurable — every rank at or above
            # The scope has all tips on one side of the split — so cases would be
            # Grading UNKNOWN markers instead of exercising the screen.
            "--gene-tree-builder": "fasttree",
            "--coalescent-mode": "off",
            "--skip-checkm": None,
            "--force": None,
        }
        # A case that passes --max-markers in `extra` overrides the default
        # Rather than duplicating it: argparse takes the last occurrence, so a
        # Duplicate would silently work while the recorded command line lied
        # About what the case asked for.
        supplied = {a for a in (str(x) for x in extra) if a.startswith("--")}
        for flag, value in defaults.items():
            if flag in omit or flag in supplied:
                continue
            args.append(flag)
            if value is not None:
                args.append(value)
        for flag, value in kwargs.items():
            flag = "--" + flag.replace("_", "-")
            if flag in omit or flag in supplied:
                continue
            if value is True:
                args.append(flag)
            elif value is not None:
                args += [flag, str(value)]
        args += [str(x) for x in extra]
        return run_markerfinder(args, expect_rc=expect_rc, timeout=timeout)

    return _mf


@pytest.fixture
def case_output(request) -> Path:
    """A clean, test-id-named output directory (the -o target)."""
    name = request.node.name.replace("/", "_").replace("::", "_")
    path = OUTPUTS / name
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    return path


@pytest.fixture
def case_workdir(request) -> Path:
    """Scratch space for files a case itself creates (configs, bad inputs)."""
    name = request.node.name.replace("/", "_").replace("::", "_")
    path = _worker_dir(WORK / "case_scratch") / name
    if path.exists():
        shutil.rmtree(path)
    path.mkdir(parents=True)
    return path


# --------------------------------------------------------------------------
# Product readers
# --------------------------------------------------------------------------
@pytest.fixture(scope="session")
def read_tsv():
    """Read a TSV written by a run into a list of dicts."""
    def _get(path: Path) -> List[Dict[str, str]]:
        lines = Path(path).read_text(encoding="utf-8").rstrip("\n").split("\n")
        header = lines[0].split("\t")
        return [dict(zip(header, line.split("\t"))) for line in lines[1:] if line]
    return _get


def _redact(value):
    """Strip machine-specific prefixes from a recorded value.

    The archived metrics ship with the release, and a path like
    ``/home/<user>/.conda/envs/markerfinder/...`` says more about the host that
    ran the suite than about the result. ``<repo>``/``<env>``/``~`` keep the
    record readable and comparable everywhere.
    """
    if isinstance(value, dict):
        return {k: _redact(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_redact(v) for v in value]
    if not isinstance(value, str):
        return value
    text = value.replace(str(REPO), "<repo>")
    home = str(Path.home())
    text = text.replace(home, "~")
    return re.sub(r"~[/\\][^\s,;)]*?[/\\]envs[/\\][A-Za-z0-9_.-]+", "<env>", text)


@pytest.fixture(scope="session")
def record_metric(results_dir: Path):
    """Append a measured number to results/metrics.json for the test report.

    The report quotes these values instead of re-typing them, so a stale
    number in prose becomes a diff in the machine-readable record.
    """
    path = results_dir / "metrics.json"

    def _record(case: str, key: str, value) -> None:
        payload = {}
        if path.exists():
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
            except json.JSONDecodeError:
                payload = {}
        payload.setdefault(case, {})[key] = _redact(value)
        path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n",
                        encoding="utf-8")

    return _record
