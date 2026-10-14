#!/usr/bin/env python3
"""Bounded, sequential, fresh-process PRIME scaling comparisons (stdlib only)."""
import argparse
import csv
import ctypes
import hashlib
import json
import os
from pathlib import Path
import platform
import signal
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def atomic_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    temporary.replace(path)


def file_tree(path):
    return {str(p.relative_to(path)): sha(p)
            for p in sorted(path.rglob("*")) if p.is_file()}


def tree_digest(files):
    return hashlib.sha256(json.dumps(files, sort_keys=True).encode()).hexdigest()


def package_snapshot(library):
    package = library / "prime"
    if not (package / "DESCRIPTION").is_file():
        raise ValueError(f"No installed prime package in {library}")
    files = file_tree(package)
    description = (package / "DESCRIPTION").read_text()
    version = next(s.split(":", 1)[1].strip()
                   for s in description.splitlines() if s.startswith("Version:"))
    return {"library": str(library), "package_version": version,
            "installed_sha256": files, "tree_sha256": tree_digest(files)}


def source_snapshot():
    names = subprocess.check_output(
        ["git", "-C", str(REPO), "ls-files", "-z"]).decode().split("\0")
    names = [s for s in names if s and (REPO / s).is_file()]
    files = {s: sha(REPO / s) for s in names}
    return {"source_repository": str(REPO), "tracked_sha256": files,
            "tracked_tree_sha256": tree_digest(files),
            "git_head": subprocess.check_output(
                ["git", "-C", str(REPO), "rev-parse", "HEAD"]).decode().strip()}


class MacTaskInfo(ctypes.Structure):
    _fields_ = [(name, ctypes.c_uint64) for name in (
        "virtual", "resident", "total_user", "total_system", "threads_user",
        "threads_system")] + [(name, ctypes.c_int32) for name in (
        "policy", "faults", "pageins", "cow_faults", "messages_sent",
        "messages_received", "syscalls_mach", "syscalls_unix", "csw",
        "threadnum", "numrunning", "priority")]


def resident_bytes(pid):
    if sys.platform == "darwin":
        info = MacTaskInfo()
        lib = ctypes.CDLL("/usr/lib/libproc.dylib", use_errno=True)
        count = lib.proc_pidinfo(pid, 4, 0, ctypes.byref(info), ctypes.sizeof(info))
        return info.resident if count == ctypes.sizeof(info) else None
    if sys.platform.startswith("linux"):
        try:
            fields = Path(f"/proc/{pid}/status").read_text().splitlines()
            line = next(s for s in fields if s.startswith("VmRSS:"))
            return int(line.split()[1]) * 1024
        except (OSError, StopIteration):
            return None
    raise RuntimeError("RSS budget enforcement supports macOS and Linux.")


def stop_process(process):
    try:
        os.killpg(process.pid, signal.SIGTERM)
        process.wait(timeout=2)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()
    except ProcessLookupError:
        process.wait()


def parse_resource_peak(stderr):
    # /usr/bin/time reports bytes on macOS and KiB on Linux.
    if sys.platform == "darwin":
        for line in stderr.splitlines():
            if "maximum resident set size" in line:
                return int(line.split()[0])
    else:
        for line in stderr.splitlines():
            if "Maximum resident set size (kbytes):" in line:
                return int(line.rsplit(":", 1)[1].strip()) * 1024
    return None


def run_one(task, output, args, libraries):
    directory = output / "runs" / task["run_id"]
    directory.mkdir(parents=True, exist_ok=False)
    snapshot = package_snapshot(libraries[task["variant"]])
    env = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                 "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        env[name] = "1"
    env["PRIME_NUM_THREADS"] = str(args.threads)
    rcommand = [args.rscript, "--vanilla", str(HERE / "worker.R"),
                str(directory), str(libraries[task["variant"]]), task["variant"],
                task["mode"], str(task["n"]), str(task["p"]), str(task["seed"]),
                str(task["neighbor_eps"]), str(args.query_block_size),
                str(args.max_selected), str(args.recall_queries),
                str(args.twin_r) if args.twin_r is not None else "NULL"]
    command = (["/usr/bin/time", "-l", "-p"] if sys.platform == "darwin"
               else ["/usr/bin/time", "-v"]) + rcommand
    invocation = dict(task, command=command, package=snapshot,
                      timeout_sec=args.timeout, rss_limit_bytes=int(args.memory_gib * 2**30),
                      threads=args.threads, blas_omp_threads=1,
                      query_block_size=args.query_block_size,
                      started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
    atomic_json(directory / "invocation.json", invocation)
    start = time.monotonic()
    peak = 0
    monitor_reads = 0
    status = "completed"
    # R writes its own PID: monitor the R process, not /usr/bin/time's small RSS.
    with (directory / "stdout.log").open("w") as stdout, \
            (directory / "stderr.log").open("w") as stderr:
        process = subprocess.Popen(command, stdout=stdout, stderr=stderr,
                                   env=env, cwd=REPO, start_new_session=True)
        try:
            while process.poll() is None:
                elapsed = time.monotonic() - start
                pid_path = directory / "worker.pid"
                try:
                    pid = int(pid_path.read_text())
                except (OSError, ValueError):
                    pid = None
                rss = resident_bytes(pid) if pid else None
                if rss is not None:
                    peak = max(peak, rss)
                    monitor_reads += 1
                if elapsed > args.timeout:
                    status = "timeout"
                    stop_process(process)
                    break
                if rss is not None and rss > args.memory_gib * 2**30:
                    status = "memory_budget_exceeded"
                    stop_process(process)
                    break
                if elapsed > 5 and monitor_reads == 0:
                    status = "rss_monitor_unavailable"
                    stop_process(process)
                    break
                time.sleep(0.10)
        except BaseException:
            stop_process(process)
            atomic_json(directory / "interrupted.json", {"reason": "runner interrupted"})
            raise
    wall = time.monotonic() - start
    peak_rusage = parse_resource_peak((directory / "stderr.log").read_text())
    row = dict(task, status=status, returncode=process.returncode, wall_sec=wall,
               peak_rss_bytes=peak_rusage, sampled_peak_rss_bytes=peak,
               rss_monitor_reads=monitor_reads, package_version=snapshot["package_version"],
               library_sha256=snapshot["tree_sha256"])
    result_path = directory / "metrics.csv"
    if result_path.exists():
        with result_path.open() as handle:
            metrics = list(csv.DictReader(handle))
        if len(metrics) != 1:
            raise RuntimeError(f"Malformed metrics: {result_path}")
        row.update(metrics[0])
    if status != "completed":
        row["status"] = status
    elif process.returncode != 0:
        row["status"] = "error"
    elif not result_path.exists():
        row["status"] = "missing_result"
    after = package_snapshot(libraries[task["variant"]])
    row["library_unchanged"] = snapshot["tree_sha256"] == after["tree_sha256"]
    if not row["library_unchanged"]:
        row["status"] = "library_changed_during_run"
    atomic_json(directory / "receipt.json", row)
    return row


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--new-lib", type=Path, required=True)
    parser.add_argument("--baseline-lib", type=Path, required=True,
                        help="Frozen baseline library supplied separately; no workstation default.")
    parser.add_argument("--baseline-build", type=Path, required=True,
                        help="Directory with the frozen baseline archive and build_manifest.json.")
    parser.add_argument("--new-build-manifest", type=Path,
                        help="Optional source/archive/installed manifest for the new build.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n", default="2000,10000,50000")
    parser.add_argument("--p", default="5,25")
    parser.add_argument("--modes", default="screen,pme")
    parser.add_argument("--eps", default="0.1,0.5")
    parser.add_argument("--variants", default="",
                        help="Optional comma-separated subset: baseline,indexed_exact,ann_eps0.1,ann_eps0.5")
    parser.add_argument("--replications", type=int, default=1)
    parser.add_argument("--seed", type=int, default=202609130)
    parser.add_argument("--threads", type=int, default=1,
                        help="PRIME candidate threads; BLAS/OMP/MKL/VECLIB remain at one.")
    parser.add_argument("--timeout", type=float, default=180)
    parser.add_argument("--total-budget", type=float, default=900,
                        help="Stop dispatching new runs after this total wall budget (seconds).")
    parser.add_argument("--memory-gib", type=float, default=4)
    parser.add_argument("--query-block-size", type=int, default=256)
    parser.add_argument("--twin-r", type=int,
                        help="Optional existing Data-Twinning inverse split ratio; default uses all queries.")
    parser.add_argument("--max-selected", type=int, default=8)
    parser.add_argument("--recall-queries", type=int, default=32)
    parser.add_argument("--rscript", default="Rscript")
    parser.add_argument("--execute", action="store_true",
                        help="Without this flag write only the immutable plan/provenance.")
    args = parser.parse_args()
    args.new_lib = args.new_lib.resolve()
    args.baseline_lib = args.baseline_lib.resolve()
    args.output = args.output.resolve()
    if args.new_lib == args.baseline_lib:
        parser.error("Baseline and new package libraries must be distinct.")
    if min(args.replications, args.threads, args.max_selected, args.recall_queries,
           args.query_block_size) < 1 or args.timeout <= 0 or args.memory_gib <= 0 or args.total_budget <= 0:
        parser.error("Counts and budgets must be positive.")
    ns = [int(s) for s in args.n.split(",")]
    ps = [int(s) for s in args.p.split(",")]
    eps = [float(s) for s in args.eps.split(",") if s]
    modes = args.modes.split(",")
    if min(ns) < 10 or min(ps) < 3 or any(e <= 0 for e in eps):
        parser.error("Require n>=10, p>=3, and positive ANN eps values.")
    if args.twin_r is not None and (args.twin_r < 2 or args.twin_r > min(ns) // 2):
        parser.error("twin-r must be an integer from 2 through floor(min(n)/2).")
    if not set(modes) <= {"screen", "pme", "pipeline", "neighbors"}:
        parser.error("Modes: screen,pme,pipeline,neighbors")
    if any(len(values) != len(set(values)) for values in (ns, ps, eps, modes)):
        parser.error("Grid values must be unique.")
    args.output.mkdir(parents=True, exist_ok=False)
    variants = {"baseline": (args.baseline_lib, 0.0),
                "indexed_exact": (args.new_lib, 0.0)}
    variants.update({f"ann_eps{e:g}": (args.new_lib, e) for e in eps})
    if args.variants:
        requested = args.variants.split(",")
        if len(requested) != len(set(requested)) or not set(requested) <= set(variants):
            parser.error("Unknown or duplicated variant.")
        variants = {name: variants[name] for name in requested}
    tasks = []
    for n in ns:
        for p in ps:
            for replication in range(1, args.replications + 1):
                for mode in modes:
                    for variant, (_, e) in variants.items():
                        if mode == "neighbors" and variant == "baseline":
                            continue  # worker uses its own exact brute-force reference.
                        tasks.append(dict(
                            run_id=f"n{n}_p{p}_r{replication}_{mode}_{variant}",
                            n=n, p=p, seed=args.seed + replication,
                            replication=replication, mode=mode, variant=variant,
                            neighbor_eps=e, twin_r=args.twin_r))
    baseline_manifest = args.baseline_build / "build_manifest.json"
    archive = args.baseline_build / "prime_0.5.3.tar.gz"
    if not baseline_manifest.is_file() or not archive.is_file():
        raise RuntimeError("Frozen baseline source archive/build manifest is missing.")
    baseline_build = json.loads(baseline_manifest.read_text())
    if sha(archive) != baseline_build["source_archive_sha256"]:
        raise RuntimeError("Frozen baseline source archive failed its recorded hash.")
    libraries = {k: v[0] for k, v in variants.items()}
    snapshots = {name: package_snapshot(lib) for name, lib in libraries.items()}
    for filename, expected in baseline_build["installed_files"].items():
        relative = Path(filename).relative_to(Path(baseline_build["library"]))
        candidate = args.baseline_lib / relative
        if sha(candidate) != expected:
            raise RuntimeError(f"Frozen baseline installed file changed: {candidate}")
    atomic_json(args.output / "baseline_build_manifest.json", baseline_build)
    new_build = None
    if args.new_build_manifest:
        new_build = json.loads(args.new_build_manifest.read_text())
        if sha(Path(new_build["archive"])) != new_build["archive_sha256"]:
            raise RuntimeError("New source archive failed its recorded hash.")
        for filename, expected in new_build["installed_sha256"].items():
            recorded = Path(filename)
            candidate = (args.new_lib / "prime" / recorded if not recorded.is_absolute()
                         else args.new_lib / recorded.relative_to(Path(new_build["installed_library"])))
            if sha(candidate) != expected:
                raise RuntimeError(f"New installed file changed: {candidate}")
        atomic_json(args.output / "new_build_manifest.json", new_build)
    manifest = dict(schema="prime-scaling-v1", tasks=tasks,
                    planned_only=not args.execute, python=sys.version,
                    platform=platform.platform(), packages=snapshots,
                    source_at_plan=source_snapshot(),
                    baseline_archive=str(archive), baseline_archive_sha256=sha(archive),
                    new_build_manifest_sha256=(sha(args.new_build_manifest)
                                               if args.new_build_manifest else None),
                    harness_sha256={p.name: sha(p) for p in HERE.iterdir()
                                    if p.is_file() and p.suffix in (".py", ".R", ".md")},
                    budgets=dict(timeout_sec=args.timeout, memory_gib=args.memory_gib,
                                 total_dispatch_budget_sec=args.total_budget),
                    threads=args.threads, blas_omp_threads=1,
                    query_block_size=args.query_block_size,
                    twin_r=args.twin_r)
    atomic_json(args.output / "manifest.json", manifest)
    print(f"Planned {len(tasks)} runs; execute={args.execute}", flush=True)
    if not args.execute:
        return
    if resident_bytes(os.getpid()) is None:
        raise RuntimeError("RSS monitor is unavailable; no benchmark launched.")
    rows = []
    dispatch_start = time.monotonic()
    for task in tasks:
        if time.monotonic() - dispatch_start >= args.total_budget:
            print("Total dispatch budget reached; remaining runs are not attempted.", flush=True)
            break
        row = run_one(task, args.output, args, libraries)
        rows.append(row)
        columns = sorted(set().union(*(set(r) for r in rows)))
        with (args.output / "raw.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns)
            writer.writeheader()
            writer.writerows(rows)
        print(task["run_id"], row["status"], f"{row['wall_sec']:.3f}s", flush=True)
    atomic_json(args.output / "dispatcher_finished.json", {
        "runs_attempted": len(rows), "runs_planned": len(tasks),
        "dispatch_complete": len(rows) == len(tasks)})


if __name__ == "__main__":
    main()
