#!/usr/bin/env python3
"""Validate lineage and compare exact/ANN outcomes without hiding incomplete runs."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys


def number(row, key):
    value = row.get(key)
    if value in (None, "", "NA"):
        return None
    return float(value)


def effects(path):
    if not path.is_file():
        return None
    rows = list(csv.DictReader(path.open()))
    return {int(r["coordinate"]): float(r["effect"]) for r in rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--exact-tolerance", type=float, default=1e-8)
    parser.add_argument("--require-all-completed", action="store_true")
    parser.add_argument("--all-query-reference", type=Path,
                        help="Saved all-query indexed-exact fixed-s PME output for query-approximation comparison.")
    args = parser.parse_args()
    output = args.output.resolve()
    manifest = json.loads((output / "manifest.json").read_text())
    receipts = {p.parent.name: json.loads(p.read_text())
                for p in (output / "runs").glob("*/receipt.json")}
    issues = []
    comparisons = []
    for task in manifest["tasks"]:
        run_id = task["run_id"]
        row = receipts.get(run_id)
        if row is None:
            comparisons.append(dict(run_id=run_id, comparison_status="not_attempted"))
            if args.require_all_completed:
                issues.append(f"Missing run: {run_id}")
            continue
        if not row.get("library_unchanged"):
            issues.append(f"Package changed during run: {run_id}")
        invocation = json.loads((output / "runs" / run_id / "invocation.json").read_text())
        if invocation["package"]["tree_sha256"] != manifest["packages"][task["variant"]]["tree_sha256"]:
            issues.append(f"Package differs from frozen plan: {run_id}")
        for key in ("n", "p", "mode", "variant", "seed"):
            if str(row[key]) != str(task[key]):
                issues.append(f"Task key mismatch {key}: {run_id}")
        recorded_twin = row.get("twin_r")
        if recorded_twin in (None, "", "NULL"):
            recorded_twin = None
        else:
            recorded_twin = int(recorded_twin)
        if recorded_twin != task.get("twin_r"):
            issues.append(f"Query-approximation setting mismatch: {run_id}")
        if row["status"] != "completed":
            comparisons.append(dict(run_id=run_id, comparison_status=row["status"]))
            if args.require_all_completed:
                issues.append(f"Incomplete run ({row['status']}): {run_id}")
            continue
        if row["returncode"] != 0 or row["rss_monitor_reads"] <= 0:
            issues.append(f"Invalid completion receipt: {run_id}")
        if row.get("actual_package_version") != invocation["package"]["package_version"]:
            issues.append(f"Loaded package version differs from installation: {run_id}")
        query_file = output / "runs" / run_id / "query_indices.txt"
        if task["mode"] != "neighbors" and query_file.exists():
            queries = [int(s) for s in query_file.read_text().splitlines()]
            if (len(queries) != int(row["query_count"]) or len(queries) != len(set(queries))
                    or not all(1 <= q <= int(task["n"]) for q in queries)):
                issues.append(f"Invalid actual query indices: {run_id}")
            if recorded_twin is None and len(queries) != int(task["n"]):
                issues.append(f"All-query default did not use every observation: {run_id}")
        if task["mode"] == "neighbors":
            ratio = number(row, "neighbor_distance_ratio_max")
            recall = number(row, "neighbor_recall_min")
            error = number(row, "neighbor_squared_distance_error_max")
            if ratio is None or ratio > 1 + task["neighbor_eps"] + 1e-9:
                issues.append(f"Neighbor distance bound failed: {run_id}")
            if task["neighbor_eps"] == 0 and (recall is None or recall < 1):
                issues.append(f"Exact neighbor recall failed: {run_id}")
            if error is None or error > 1e-8:
                issues.append(f"Returned squared distances failed: {run_id}")
            comparisons.append(dict(run_id=run_id, comparison_status="neighbor_reference_checked",
                                    neighbor_recall_min=recall, neighbor_distance_ratio_max=ratio))
            continue
        pair_key = {k: task[k] for k in ("n", "p", "replication", "mode", "seed")}
        references = [t for t in manifest["tasks"]
                      if all(t[k] == v for k, v in pair_key.items())
                      and t["variant"] in ("baseline", "indexed_exact")]
        for reference_task in references:
            if reference_task["run_id"] == run_id:
                continue
            reference = receipts.get(reference_task["run_id"])
            if reference is None or reference["status"] != "completed":
                comparisons.append(dict(run_id=run_id, reference=reference_task["variant"],
                                        comparison_status="reference_incomplete"))
                continue
            if row.get("input_md5") != reference.get("input_md5"):
                issues.append(f"Paired inputs differ: {run_id}")
            reference_query = output / "runs" / reference_task["run_id"] / "query_indices.txt"
            if query_file.exists() and reference_query.exists() and query_file.read_bytes() != reference_query.read_bytes():
                issues.append(f"Paired query locations differ: {run_id}")
            same_candidates = row.get("candidates") == reference.get("candidates")
            same_selected = row.get("selected") == reference.get("selected")
            current_effects = effects(output / "runs" / run_id / "effects.csv")
            reference_effects = effects(output / "runs" / reference_task["run_id"] / "effects.csv")
            half_l1 = None
            if current_effects is not None and reference_effects is not None:
                if set(current_effects) != set(reference_effects):
                    issues.append(f"Allocation dimensions differ: {run_id}")
                else:
                    half_l1 = sum(abs(current_effects[k] - reference_effects[k])
                                  for k in current_effects) / 2
                for label, values in ((run_id, current_effects), (reference_task["run_id"], reference_effects)):
                    if not all(math.isfinite(v) and v >= -args.exact_tolerance for v in values.values()):
                        issues.append(f"Invalid allocation: {label}")
                    if abs(sum(values.values()) - 1) > args.exact_tolerance and any(values.values()):
                        issues.append(f"Allocation not normalized: {label}")
            exact_pair = task["neighbor_eps"] == 0 and reference_task["neighbor_eps"] == 0
            if exact_pair and (not same_candidates or not same_selected or
                               (half_l1 is not None and half_l1 > args.exact_tolerance)):
                issues.append(f"Exact implementation parity failed: {run_id}")
            current_time = number(row, "operation_sec")
            reference_time = number(reference, "operation_sec")
            screen_time = number(row, "screen_sec")
            reference_screen_time = number(reference, "screen_sec")
            pme_time = number(row, "pme_sec")
            reference_pme_time = number(reference, "pme_sec")
            comparisons.append(dict(
                run_id=run_id, reference=reference_task["variant"], comparison_status="paired_completed",
                candidates_identical=same_candidates, selected_identical=same_selected,
                allocation_half_l1=half_l1,
                operation_speedup=(reference_time / current_time if current_time and reference_time else None),
                screen_speedup=(reference_screen_time / screen_time
                                if screen_time and reference_screen_time else None),
                pme_speedup=(reference_pme_time / pme_time
                             if pme_time and reference_pme_time else None),
                peak_rss_bytes=row.get("peak_rss_bytes"),
                reference_peak_rss_bytes=reference.get("peak_rss_bytes")))
    if args.all_query_reference:
        reference_output = args.all_query_reference.resolve()
        reference_rows = [json.loads(p.read_text())
                          for p in (reference_output / "runs").glob("*/receipt.json")]
        for task in manifest["tasks"]:
            row = receipts.get(task["run_id"])
            if row is None or row["status"] != "completed" or task["mode"] not in ("pipeline", "pme"):
                continue
            matched = [r for r in reference_rows if r["status"] == "completed"
                       and r["variant"] == "indexed_exact" and r["mode"] == "pme"
                       and r.get("twin_r") in (None, "", "NULL")
                       and all(str(r[k]) == str(task[k]) for k in ("n", "p", "seed", "replication"))]
            if len(matched) != 1:
                issues.append(f"Expected one completed all-query exact PME reference: {task['run_id']}")
                continue
            reference = matched[0]
            if row.get("input_md5") != reference.get("input_md5"):
                issues.append(f"Inputs differ from all-query reference: {task['run_id']}")
            same_selected = row.get("selected") == reference.get("selected")
            current_effects = effects(output / "runs" / task["run_id"] / "effects.csv")
            reference_effects = effects(reference_output / "runs" / reference["run_id"] / "effects.csv")
            distance = None
            if current_effects is not None and reference_effects is not None and set(current_effects) == set(reference_effects):
                distance = sum(abs(current_effects[k] - reference_effects[k]) for k in current_effects) / 2
            comparisons.append(dict(
                run_id=task["run_id"], reference="saved_all_query_indexed_exact_fixed_s_PME",
                comparison_status="paired_completed", selected_identical=same_selected,
                allocation_half_l1=distance, actual_query_count=row.get("query_count"),
                reference_query_count=reference.get("query_count"),
                approximation_sources=("query only" if task["neighbor_eps"] == 0 else "query and ANN"),
                reference_output=str(reference_output)))
    columns = sorted(set().union(*(set(r) for r in comparisons))) if comparisons else ["run_id"]
    with (output / "comparisons.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(comparisons)
    report = dict(schema="prime-scaling-validation-v1", planned=len(manifest["tasks"]),
                  validator_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  attempted=len(receipts), completed=sum(r["status"] == "completed" for r in receipts.values()),
                  exact_tolerance=args.exact_tolerance, issues=issues,
                  evidence_valid=not issues and bool(receipts),
                  completed_evidence_exists=any(r["status"] == "completed" for r in receipts.values()),
                  all_runs_completed=len(receipts) == len(manifest["tasks"]) and
                  all(r["status"] == "completed" for r in receipts.values()))
    (output / "validation.json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))
    return 1 if issues else 0


if __name__ == "__main__":
    sys.exit(main())
