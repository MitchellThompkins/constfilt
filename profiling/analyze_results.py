#!/usr/bin/env -S uv run
# /// script
# dependencies = ["matplotlib"]
# ///
"""Analyze profiling CSV results, print summary tables, and generate plots.

Handles two CSV types (auto-detected from header):
  compile-time:  filter_type,order,method,compile_time_sec,max_rss_kb,exit_code
  runtime:       library,filter_type,order,method,ns_per_sample,msa_per_s,dc_gain
"""

import csv
import os
import sys
from collections import defaultdict


def read_metadata(csv_path):
    family, compiler, os_name = None, None, None
    with open(csv_path) as f:
        for line in f:
            if line.startswith("# family:"):
                family = line[len("# family:"):].strip()
            elif line.startswith("# compiler:"):
                compiler = line[len("# compiler:"):].strip()
            elif line.startswith("# os:"):
                os_name = line[len("# os:"):].strip()
            elif not line.startswith("#"):
                break
    return family, compiler, os_name


def detect_mode(csv_path):
    with open(csv_path) as f:
        for line in f:
            if not line.startswith("#"):
                if "ns_per_sample" in line:
                    return "runtime"
                if "max_b_err" in line:
                    return "accuracy"
                return "compiletime"
    return "compiletime"


# ── Compile-time ──────────────────────────────────────────────────────────────

def load_compiletime(csv_path):
    success = defaultdict(list)  # (filter_type, method, order) -> [sec, ...]
    failed  = defaultdict(list)
    memory  = defaultdict(list)
    orders, categories = set(), set()

    with open(csv_path) as f:
        lines = [l for l in f if not l.startswith("#")]

    for row in csv.DictReader(lines):
        order  = int(row["order"])
        ftype  = row["filter_type"]
        method = row["method"]
        cat    = f"{ftype}_{method}"
        t      = float(row["compile_time_sec"])
        orders.add(order)
        categories.add(cat)
        key = (cat, order)
        if int(row["exit_code"]) == 0:
            success[key].append(t)
            memory[key].append(float(row["max_rss_kb"]))
        elif int(row["exit_code"]) != 124:
            failed[key].append(t)

    return success, failed, memory, sorted(orders), sorted(categories)


def print_compiletime_table(success, failed, orders, categories):
    header = f"{'category':<30}" + "".join(f"{s:>8}" for s in orders)
    print(header)
    print("-" * len(header))
    for cat in categories:
        row = f"{cat:<30}"
        for o in orders:
            sv = success.get((cat, o), [])
            fv = failed.get((cat, o), [])
            if sv:
                row += f"{sum(sv)/len(sv):>8.2f}"
            elif fv:
                row += f"{'F'+f'{sum(fv)/len(fv):.1f}':>8}"
            else:
                row += f"{'---':>8}"
        print(row)
    print()
    print("Values: mean compile time in seconds. F<t> = compiler limit hit.")
    print(f"OK: {sum(len(v) for v in success.values())}  "
          f"Failed: {sum(len(v) for v in failed.values())}")


def plot_compiletime(success, failed, orders, categories, csv_path, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import itertools

    fig, ax = plt.subplots(figsize=(12, 7))
    colors = {c: col["color"] for c, col in
              zip(categories, itertools.cycle(plt.rcParams["axes.prop_cycle"]))}

    for cat in categories:
        xs = [o for o in orders if success.get((cat, o))]
        ys = [sum(success[(cat, o)]) / len(success[(cat, o)]) for o in xs]
        if xs:
            ax.plot(xs, ys, marker="o", label=cat, color=colors[cat])
        fxs = [o for o in orders if failed.get((cat, o))]
        fys = [sum(failed[(cat, o)]) / len(failed[(cat, o)]) for o in fxs]
        if fxs:
            ax.plot(fxs, fys, marker="x", markersize=10, linestyle="none",
                    color=colors[cat])

    ax.set_xlabel("Filter Order")
    ax.set_ylabel("Mean Compile Time (s)")
    ax.set_title(f"constfilt Compile-Time Cost by Filter Order\n{label}")
    ax.set_xticks(orders)
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    fig.tight_layout()

    png_path = os.path.splitext(csv_path)[0] + ".png"
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    print(f"Plot saved: {png_path}")


# ── Runtime ───────────────────────────────────────────────────────────────────

def load_runtime(csv_path):
    data   = defaultdict(list)   # (library, filter_type, method, order) -> [ns, ...]
    dc     = defaultdict(list)   # same key -> [dc_gain, ...]
    orders, groups = set(), set()

    with open(csv_path) as f:
        lines = [l for l in f if not l.startswith("#")]

    for row in csv.DictReader(lines):
        order  = int(row["order"])
        lib    = row["library"]
        ftype  = row["filter_type"]
        method = row["method"]
        grp    = f"{lib}_{ftype}_{method}"
        key    = (grp, order)
        orders.add(order)
        groups.add(grp)
        data[key].append(float(row["ns_per_sample"]))
        dc[key].append(float(row["dc_gain"]))

    return data, dc, sorted(orders), sorted(groups)


def print_runtime_table(data, dc, orders, groups):
    # Determine which libraries/groups support which methods.
    # iir1 and KFR don't have separate ZOH/MatchedZ implementations;
    # their method label is "runtime" / "runtime+simd" respectively.
    # Show "N/A" in cells where a library simply has no data for a method.
    constfilt_methods = {"tustin", "zoh", "matchedz"}
    all_methods = sorted({grp.rsplit("_", 1)[-1] for grp in groups})

    # Throughput table
    header = f"{'group':<42}" + "".join(f"{o:>8}" for o in orders)
    print("Throughput (ns/sample)  [N/A = method not supported by this library]:")
    print(header)
    print("-" * len(header))
    for grp in groups:
        row = f"{grp:<42}"
        for o in orders:
            vals = data.get((grp, o), [])
            if vals:
                row += f"{sum(vals)/len(vals):>8.2f}"
            else:
                # Distinguish "no data for this order" from "method not applicable"
                lib = grp.split("_")[0]
                method = grp.rsplit("_", 1)[-1]
                if lib != "constfilt" and method in constfilt_methods:
                    row += f"{'N/A':>8}"
                else:
                    row += f"{'---':>8}"
        print(row)
    print()

    # DC gain table
    print("DC gain (should be ~1.0 for lowpass):")
    print(header)
    print("-" * len(header))
    for grp in groups:
        row = f"{grp:<42}"
        for o in orders:
            vals = dc.get((grp, o), [])
            row += f"{sum(vals)/len(vals):>8.5f}" if vals else f"{'---':>8}"
        print(row)
    print()


def plot_runtime(data, orders, groups, csv_path, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import itertools

    fig, ax = plt.subplots(figsize=(12, 7))
    colors = {g: c["color"] for g, c in
              zip(groups, itertools.cycle(plt.rcParams["axes.prop_cycle"]))}

    for grp in groups:
        xs = [o for o in orders if data.get((grp, o))]
        ys = [sum(data[(grp, o)]) / len(data[(grp, o)]) for o in xs]
        if xs:
            # Mark KFR visually to signal different measurement methodology
            ls = "--" if "kfr" in grp else "-"
            display = grp + (" (batch=256)" if "kfr" in grp else "")
            ax.plot(xs, ys, marker="o", linestyle=ls, label=display,
                    color=colors[grp])

    ax.set_xlabel("Filter Order")
    ax.set_ylabel("ns / sample (lower = faster)")
    ax.set_title(
        f"Filter Runtime Throughput\n{label}\n"
        "Note: KFR uses SIMD batch processing (dashed); others are single-sample"
    )
    ax.set_xticks(orders)
    ax.grid(True, alpha=0.3)
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    fig.tight_layout()

    png_path = os.path.splitext(csv_path)[0] + ".png"
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    print(f"Plot saved: {png_path}")


# ── Accuracy ─────────────────────────────────────────────────────────────────

def load_accuracy(csv_path):
    data   = defaultdict(dict)   # group -> {order: max_step_err}
    orders, groups = set(), set()

    with open(csv_path) as f:
        lines = [l for l in f if not l.startswith("#")]

    for row in csv.DictReader(lines):
        order  = int(row["order"])
        lib    = row["library"]
        ftype  = row["filter_type"]
        method = row["method"]
        grp    = f"{lib}_{ftype}_{method}"
        orders.add(order)
        groups.add(grp)
        step_str = row["max_step_err"]
        step_val = float(step_str) if step_str not in ("n/a", "") else None
        b_str    = row["max_b_err"]
        b_val    = float(b_str)   if b_str    not in ("n/a", "") else None
        a_str    = row["max_a_err"]
        a_val    = float(a_str)   if a_str    not in ("n/a", "") else None
        data[grp][order] = {"b": b_val, "a": a_val, "step": step_val}

    return data, sorted(orders), sorted(groups)


ACCURACY_THRESHOLD = 1e-6  # step errors above this are flagged


def _fmt_err(val):
    if val > ACCURACY_THRESHOLD:
        return f"!{val:.0e}"[1:]    # e.g. "!2e-5" -> flag it
    return f"{val:.0e}"


def print_accuracy_table(data, orders, groups):
    print(f"Max step-response error vs Octave  (! = exceeds {ACCURACY_THRESHOLD:.0e}):\n")
    header = f"{'group':<28}" + "".join(f"{o:>9}" for o in orders)
    print(header)
    print("-" * len(header))
    for grp in groups:
        row_str = f"{grp:<28}"
        for o in orders:
            entry = data[grp].get(o)
            if entry is None:
                row_str += f"{'---':>9}"
            elif entry["step"] is None:
                row_str += f"{'N/A':>9}"
            else:
                v = entry["step"]
                cell = f"!{v:.1e}" if v > ACCURACY_THRESHOLD else f"{v:.1e}"
                row_str += f"{cell:>9}"
        print(row_str)
    print()
    print(f"Threshold: {ACCURACY_THRESHOLD:.0e}  (matches constfilt test suite tolerance)")


def plot_accuracy(data, orders, groups, csv_path, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import itertools

    fig, ax = plt.subplots(figsize=(12, 7))
    colors = {g: c["color"] for g, c in
              zip(groups, itertools.cycle(plt.rcParams["axes.prop_cycle"]))}

    for grp in groups:
        xs = [o for o in sorted(data[grp].keys())
              if data[grp][o]["step"] is not None]
        ys = [data[grp][o]["step"] for o in xs]
        if xs:
            ax.semilogy(xs, ys, marker="o", label=grp, color=colors[grp])

    ax.axhline(ACCURACY_THRESHOLD, color="red", linestyle="--",
               label=f"threshold ({ACCURACY_THRESHOLD:.0e})")
    ax.set_xlabel("Filter Order")
    ax.set_ylabel("Max |step error| vs Octave")
    ax.set_title(f"constfilt Accuracy by Filter Order and Method\n{label}")
    ax.set_xticks(orders)
    ax.grid(True, alpha=0.3, which="both")
    ax.legend(bbox_to_anchor=(1.05, 1), loc="upper left", fontsize=8)
    fig.tight_layout()

    png_path = os.path.splitext(csv_path)[0] + ".png"
    fig.savefig(png_path, dpi=150)
    plt.close(fig)
    print(f"Plot saved: {png_path}")


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    if len(sys.argv) < 2:
        print("Usage: analyze_results.py <results.csv>")
        sys.exit(1)

    csv_path = sys.argv[1]
    mode = detect_mode(csv_path)
    family, compiler, os_name = read_metadata(csv_path)

    basename = os.path.splitext(os.path.basename(csv_path))[0]
    label = basename.replace("_", " ", 1)
    if os_name:
        label += f" ({os_name})"

    print(f"Mode:     {mode}")
    print(f"Compiler: {compiler or 'unknown'}")
    print()

    if mode == "compiletime":
        success, failed, memory, orders, categories = load_compiletime(csv_path)
        print_compiletime_table(success, failed, orders, categories)
        try:
            plot_compiletime(success, failed, orders, categories, csv_path, label)
        except ImportError:
            print("(matplotlib not available, skipping plot)")
    elif mode == "accuracy":
        data, orders, groups = load_accuracy(csv_path)
        print_accuracy_table(data, orders, groups)
        try:
            plot_accuracy(data, orders, groups, csv_path, label)
        except ImportError:
            print("(matplotlib not available, skipping plot)")
    else:
        data, dc, orders, groups = load_runtime(csv_path)
        print_runtime_table(data, dc, orders, groups)
        try:
            plot_runtime(data, orders, groups, csv_path, label)
        except ImportError:
            print("(matplotlib not available, skipping plot)")


if __name__ == "__main__":
    main()
