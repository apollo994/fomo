#!/usr/bin/env python3
"""Audit requested vs. actually-used resources in a Nextflow execution trace.

Usage:
    python3 13_resource_rightsizing_audit.py <path/to/execution_trace.txt> [...]

Companion to plans/13_resource_rightsizing.md — this is the script that produced the
numbers in that plan, and the one to re-run for its verification protocol.

Interpretation notes:
  * `peak_rss` for tasks shorter than ~5 s is under-sampled by Nextflow's trace poller.
    The RELIABLE section at the end restricts to realtime >= 5 s; trust that one when
    deciding a new tier.
  * `CACHED` rows carry the metrics of their ORIGINAL execution, so a resumed run's trace
    describes several runs at once. Validate tier changes against a fresh run
    (new -w workdir), not a resumed one.
  * GB.h / cpu.h are reserved-vs-used integrals over each task's realtime. They rank where
    waste actually costs something, as opposed to where the ratio merely looks alarming.
"""

import csv
import re
import sys
from collections import defaultdict

_MEM = {"B": 1 / 1024**3, "KB": 1 / 1024**2, "MB": 1 / 1024, "GB": 1.0, "TB": 1024.0}
_DUR = {"ms": 0.001, "s": 1.0, "m": 60.0, "h": 3600.0, "d": 86400.0}


def gb(value):
    """'6 GB' / '5.4 MB' -> float GiB. Missing or '-' -> 0.0."""
    m = re.match(r"([\d.]+)\s*([KMGT]?B)", (value or "").strip())
    return float(m.group(1)) * _MEM[m.group(2)] if m else 0.0


def secs(value):
    """'1h 2m 3s' / '32.2s' / '10m' -> float seconds. Missing or '-' -> 0.0."""
    return sum(
        float(n) * _DUR[u] for n, u in re.findall(r"([\d.]+)\s*(ms|s|m|h|d)", value or "")
    )


def cores(value):
    """'68.3%' -> 0.683 cores actually consumed (Nextflow reports cumulative %CPU)."""
    return float((value or "0").rstrip("%") or 0) / 100.0


def audit(path):
    with open(path) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        print(f"{path}: empty trace")
        return

    status = defaultdict(int)
    attempts = defaultdict(int)
    for r in rows:
        status[r["status"]] += 1
        attempts[r["attempt"]] += 1

    print(f"\n{'=' * 118}\n{path}\n{'=' * 118}")
    print(f"tasks: {len(rows)}   status: {dict(status)}   attempts: {dict(attempts)}")
    if any(a != "1" for a in attempts):
        print("  !! retries present — a tier may be cut too far; check the `exit` column "
              "(137=OOM, 143=timeout)")

    res_gbh = used_gbh = res_cpuh = used_cpuh = 0.0
    wall = queue = 0.0
    per = defaultdict(lambda: defaultdict(float))
    peaks = defaultdict(list)

    for r in rows:
        rt = secs(r["realtime"])
        h = rt / 3600.0
        req_m, pk_m = gb(r["memory"]), gb(r["peak_rss"])
        req_c, use_c = float(r["cpus"]), cores(r["%cpu"])
        res_gbh += req_m * h
        used_gbh += pk_m * h
        res_cpuh += req_c * h
        used_cpuh += use_c * h
        wall += rt
        queue += max(0.0, secs(r["duration"]) - rt)

        k = r["process"].replace("FOMO:", "")
        p = per[k]
        p["n"] += 1
        p["res_gbh"] += req_m * h
        p["used_gbh"] += pk_m * h
        p["res_cpuh"] += req_c * h
        p["used_cpuh"] += use_c * h
        p["req_mem"] = req_m
        p["req_cpu"] = req_c
        p["req_time"] = secs(r["time"])
        p["max_peak"] = max(p["max_peak"], pk_m)
        p["max_cores"] = max(p["max_cores"], use_c)
        p["max_real"] = max(p["max_real"], rt)
        if rt >= 5.0:
            peaks[k].append((pk_m, use_c, rt, r["memory"], r["cpus"], r["time"]))

    def ratio(a, b):
        return f"{a / b:.1f}x" if b else "  inf"

    print(
        f"\nreserved memory {res_gbh:8.2f} GB.h   used {used_gbh:6.2f} GB.h   "
        f"overestimate {ratio(res_gbh, used_gbh)}"
    )
    print(
        f"reserved cpu    {res_cpuh:8.2f} cpu.h  used {used_cpuh:6.2f} cpu.h  "
        f"overestimate {ratio(res_cpuh, used_cpuh)}"
    )
    print(f"actual compute  {wall / 3600:8.2f} h      queue+poll overhead {queue / 3600:6.2f} h")

    print(f"\n{'-- per process, ranked by GB.h wasted ' + '-' * 80}")
    hdr = (f"{'process':46s} {'n':>4} {'req mem':>8} {'maxpeak':>8} {'mem x':>7} "
           f"{'cpus':>5} {'maxcore':>7} {'cpu x':>6} {'req t':>7} {'max t':>7} "
           f"{'GBh wast':>8}")
    print(hdr)
    for k, p in sorted(per.items(), key=lambda kv: -(kv[1]["res_gbh"] - kv[1]["used_gbh"])):
        print(
            f"{k[:46]:46s} {int(p['n']):4d} {p['req_mem']:7.1f}G {p['max_peak']:7.2f}G "
            f"{ratio(p['req_mem'], p['max_peak']):>7} {p['req_cpu']:5.0f} "
            f"{p['max_cores']:7.2f} {ratio(p['req_cpu'], p['max_cores']):>6} "
            f"{p['req_time'] / 60:6.0f}m {p['max_real']:6.0f}s "
            f"{p['res_gbh'] - p['used_gbh']:8.2f}"
        )
    print(
        f"{'TOTAL':46s} {len(rows):4d} {'':8} {'':8} {'':7} {'':5} {'':7} {'':6} "
        f"{'':7} {'':7} {res_gbh - used_gbh:8.2f}"
    )

    reliable = sum(len(v) for v in peaks.values())
    print(
        f"\n{'-- RELIABLE subset: realtime >= 5s ' + '-' * 83}\n"
        f"{reliable}/{len(rows)} tasks measured long enough to trust peak_rss. "
        f"Base tier decisions on these."
    )
    print(f"{'process':46s} {'n':>4} {'request':>18} {'maxpeak':>8} {'maxcore':>8} {'maxreal':>8}")
    for k, v in sorted(peaks.items(), key=lambda kv: -len(kv[1])):
        req = f"{v[0][3]}/{v[0][4]}cpu/{v[0][5]}"
        print(
            f"{k[:46]:46s} {len(v):4d} {req:>18} {max(x[0] for x in v):7.2f}G "
            f"{max(x[1] for x in v):8.2f} {max(x[2] for x in v):7.0f}s"
        )


if __name__ == "__main__":
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    for arg in sys.argv[1:]:
        audit(arg)
