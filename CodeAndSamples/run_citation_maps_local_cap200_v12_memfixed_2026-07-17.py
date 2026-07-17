#!/usr/bin/env python3
"""
Multi-Hop Incoming Citation Network Builder — LOCAL script version
Predatory Publishing Citation Mapping Study (Peterson et al.)

Runs locally (not Colab) so it isn't subject to session/idle limits.
Fetches citation data from OpenAlex, expanding IN (incoming citation)
networks up to IN_DEPTH hops for two seed DOI sets (SetA = predatory,
SetB = INANE/legitimate).

Designed to run unattended overnight:
  - Uses a thread pool to issue several OpenAlex requests concurrently
    (the bottleneck is network latency, not CPU, so this is the main
    lever for speed).
  - Prints progress continuously (targets processed, elapsed time,
    calls/sec, running node/edge counts) to both the console and a
    log file, so you can check on it at any time with:
        tail -f run_log_SetA.txt
  - Checkpoints partial results to CSV every CHECKPOINT_EVERY targets,
    so an interruption doesn't lose more than a few minutes of progress.

USAGE
-----
    python run_citation_maps_local.py --set A         # SetA only
    python run_citation_maps_local.py --set B         # SetB only
    python run_citation_maps_local.py --set both       # both, sequentially
    python run_citation_maps_local.py --set both --depth 2   # override depth

Leave it running overnight with, e.g. on macOS/Linux:
    nohup python run_citation_maps_local.py --set both > overnight.log 2>&1 &
    # then check progress any time with:
    tail -f overnight.log
"""

import argparse
import collections
import csv
import heapq
import json
import os
import random
import re
import shutil
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from threading import Lock

import matplotlib
matplotlib.use("Agg")  # no display needed — just save PNGs
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
import networkx as nx
import pandas as pd
import requests

# ============================================================
# USER PARAMETERS
# ============================================================

IN_DEPTH_DEFAULT = 3            # 1 = direct citers; 2 = citers-of-citers; 3 = one more
MAX_CITERS_PER_NODE = 200       # raised from 20 — the previous cap was hit by EVERY top-50 hub
                                 # node in both sets, meaning true hub connectivity was truncated
                                 # identically in both groups, likely compressing any real
                                 # difference between them. 200 is OpenAlex's max per_page in a
                                 # single request (no pagination needed).
NODES_HARD_CAP_PER_SET = 200000  # raised from 60,000 to give headroom for the cap-200 run

CONTACT_EMAIL = "gpeterson@nccu.edu"

OPENALEX_API_KEY = "4SJDaluvTK6Bzl0HZGLqrU"

MAX_WORKERS = 4                 # lowered from 8 for tonight's overnight run — tonight's SetInane
                                 # attempt hit sustained 429s even at lower volume; fewer parallel
                                 # requests reduces the odds of retriggering that. Can raise back up
                                 # via --workers if it turns out to be too conservative.
MAX_RETRIES = 8                 # raised from 4 — tonight showed sustained (not transient) 429s

CHECKPOINT_EVERY = 25           # save partial CSVs every N targets processed
PROGRESS_EVERY = 10             # print a progress line every N targets processed

COLOR_BY_DEPTH = True
MAX_NODES_TO_DRAW = 1200        # cap for plotting only (doesn't affect data)

OUT_BASE_DIR_IN = "citation_maps_output_local_cap200"

LABEL_A = "SetA"   # Predatory / Cabell's sample (PP)
LABEL_B = "SetB"   # Legitimate / INANE sample

# --- Full N=60 seed DOIs, from PPCitationSamples_PY_completed.csv ---
DOIS_SET_A = [
    "10.46882/AJNM/1236",
    "10.12691/ajnr-13-4-6",
    "10.11648/j.ajns.20251404.12",
    "10.25159/2520-5293/17183",
    "10.26420/annnursrespract.2023.1058",
    "10.20431/2455-4324.1101001",
    "10.17352/2581-4265.000058",
    "10.52711/2349-2996.2025.00040",
    "10.24203/ncxxjr23",
    "10.7537/marsbnj110425.01",
    "10.5430/cns.v14n1p6",
    "http://dx.doi.org/10.31031/cojnh.2025.09.000712",
    "10.54026/CFJN/1001",
    "10.36648/2049-5471.22.1.57",
    "10.36349/easjnm.2026.v08i01.003",
    "10.28933/gjn-2021-11-0505",
    "10.36648/1791-809X.18.12.1205",
    "10.53555/hsn.v12i1.6549",
    "10.47310/iajapn.2025.v06i01.001",
    "http://doi.org/10.23937/2469-5823/1510207",
    "10.23937/2469-5823/1510203",
    "10.53555/hsn.v11i1.2449",
    "10.52711/2454-2652.2025.00046",
    "10.15640/ijn.v12p1",
    "10.31690/ijnmi.2025.v010i04.007",
    "10.37745/ijnmh.15/vol12n1132",
    "10.46610/IJRMSN.2023.v04i02.001",
    "http://dx.doi.org/10.33552/IJNC.2025.05.000619",
    "10.37421/2573-0347.2024.9.365",
    "10.4172/2471-9846.1000207",
    "10.46610/JCSHN.2023.v05i02.001",
    "10.7176/JHMN/119-08",
    "10.36266/JHNR/101",
    "10.36846/2472-1654-10.1.64",
    "10.35248/2471-8505-10.1.01",
    "10.46610/JMWHGN.2026.v08i01.001",
    "10.46610/JNVI.2025.v07i03.003",
    "10.46610/JNVI.2026.v08i01.001",
    "10.37421/2167-1168.2025.14.686",
    "http://doi.org/10.36959/545/430",
    "10.29011/2577-1450.1000202",
    "10.47485/3065-7636.1043",
    "10.47485/3065-7636.1043",
    "10.5430/jnep.v15n6p1",
    "10.35248/2573-4598.23.9.248",
    "10.24966/PPN-5681/100061",
    "10.14744/phd.2025.92603",
    "10.14303/2315-568X.2025.66",
    "10.18689/mjn-1000134",
    "10.23880/nhij-16000109",
    "10.13189/nh.2025.100101",
    "10.22259/2639-1783.0801004",
    "10.4236/ojn.2026.162009",
    "10.36648/1479-1064.32.4.20",
    "10.3389/fpsyg.2018.00348",
    "10.36348/sijcms.2026.v09i01.005",
    "http://dx.doi.org/10.15226/2471-6529/7/1/00151",
    "10.14445/24547484/IJNHS-V11I3P101",
    "http://dx.doi.org/10.2174/0118744346378124250529095629",
    "http://dx.doi.org/10.51521/WJNP.2025.11101",
]

DOIS_SET_B = [
    "10.1097/TME.0000000000000370",
    "10.1097/ANC.0000000000001344",
    "10.4037/ajcc2026265",
    "10.51256/ANJ022628",
    "10.18776/95jvxc48",
    "10.1016/j.auec.2025.06.010",
    "10.2478/ajon-2025-0012",
    "10.1016/j.aucc.2025.101484",
    "10.1097/NCC.0000000000001549",
    "10.15452/cejnm.2025.16.0020",
    "10.1188/26.CJON.41-46",
    "10.1097/NUR.0000000000000940",
    "10.1080/24694193.2025.2524671",
    "10.4037/ccn2026884",
    "10.1016/j.ccell.2011.06.002",
    "10.1136/ebnurs-2025-104480",
    "10.1097/SGA.0000000000000919",
    "10.1097/HMR.0000000000000458",
    "10.1177/19375867251406198",
    "10.1016/j.ijans.2025.100917",
    "10.1097/XEB.0000000000000486",
    "10.1016/j.ijotn.2026.101257",
    "10.1111/jspn.70010",
    "10.1111/jan.70491",
    "10.1891/JDNP-2021-0023",
    "10.1177/0898010111412189",
    "10.1097/JNN.0000000000000854",
    "10.1097/nna.0000000000000824",
    "10.1097/NCQ.0000000000000892",
    "10.13178/jnparr.2025.1501.1503",
    "10.1097/jnr.0000000000000311",
    "10.1016/j.jogn.2026.01.005",
    "10.1016/j.pedhc.2025.08.009",
    "10.1177/23320249261417953",
    "10.1097/JPN.0000000000000960",
    "10.1177/10598405251409622",
    "10.1097/JXX.0000000000001238",
    "10.2309/1557-1289-30.4.42",
    "10.1097/JNC.0000000000000564",
    "10.1177/10436596251370370",
    "10.1097/JTN.0000000000000882",
    "10.7454/jki.v28i3.1270",
    "10.1097/NMC.0000000000001152",
    "10.62116/MSJ.2025.34.6.266",
    "10.7748/nr.2025.e1967",
    "10.1111/nhs.70303",
    "10.1111/nicc.70356",
    "10.1016/j.outlook.2025.102523",
    "10.7748/ns2002.10.17.5.33.c3282",
    "10.62116/PNJ.2025.51.6.265",
    "10.1177/15269248251383949",
    "10.55048/jpns194",
    "10.1016/j.wcn.2025.05.002",
    "10.1111/wvn.70083",
    "10.7748/mhp.2025.e1755",
    "10.1111/jorc.70051",
    "10.1016/j.ijnsa.2025.100443",
    "10.24198/jkp.v13i3.2748",
    "10.1515/ijnes-2025-0051",
    "10.1177/10848223241308483",
]

# SetA (predatory) is left at the full 60 DOIs — natural attrition (~10 unresolvable
# DOIs) determines which ~50 actually resolve, same as before; we don't know in advance
# which ones will fail, so we don't pre-select. SetB (INANE) is truncated to the first 50
# to roughly match SetA's typical resolved count for a comparable sample size.
DOIS_SET_B = DOIS_SET_B[:50]

# ============================================================
# OpenAlex helpers
# ============================================================

def oa_get(url, params=None):
    """GET with OpenAlex polite-pool identification and exponential backoff retry.
    Handles SUSTAINED rate-limiting (many consecutive 429s), not just transient ones —
    backoff grows to several minutes with jitter, and a real exception is always raised
    on exhaustion (never None), since a bare 'raise None' crashes with a confusing,
    unrelated-looking error — exactly what killed tonight's SetInane run."""
    params = dict(params or {})
    params.setdefault("mailto", CONTACT_EMAIL)
    params.setdefault("api_key", OPENALEX_API_KEY)
    last_err = RuntimeError(f"oa_get exhausted all {MAX_RETRIES} retries for {url} "
                             f"(likely sustained rate-limiting)")
    for attempt in range(MAX_RETRIES):
        try:
            r = requests.get(url, params=params, timeout=45)
            if r.status_code == 429:
                wait = min((2 ** attempt) * 5 + random.uniform(0, 3), 300)  # cap at 5 min, with jitter
                log(f"  [429 rate limited] backing off {wait:.0f}s (attempt {attempt+1}/{MAX_RETRIES})...")
                time.sleep(wait)
                continue
            r.raise_for_status()
            return r.json()
        except requests.exceptions.HTTPError as e:
            status = e.response.status_code if e.response is not None else None
            if status in (404, 400, 403):
                # Permanent failure — this DOI/request will never succeed no matter how
                # many times we retry it. Fail immediately instead of burning minutes.
                log(f"  [{status} — not retrying, permanent failure] {e}")
                raise
            last_err = e
            wait = min((2 ** attempt) * 2 + random.uniform(0, 2), 120)
            log(f"  [HTTP error {status}] retrying in {wait:.0f}s (attempt {attempt+1}/{MAX_RETRIES})...")
            time.sleep(wait)
        except Exception as e:
            last_err = e
            wait = min((2 ** attempt) * 2 + random.uniform(0, 2), 120)
            log(f"  [non-429 error: {type(e).__name__}: {e}] retrying in {wait:.0f}s "
                f"(attempt {attempt+1}/{MAX_RETRIES})...")
            time.sleep(wait)
    raise last_err


def work_from_doi(doi):
    url = f"https://api.openalex.org/works/https://doi.org/{doi.strip().lower()}"
    return oa_get(url)


def get_citers(openalex_id, per_page=25):
    url = "https://api.openalex.org/works"
    return oa_get(url, params={"filter": f"cites:{openalex_id}", "per_page": per_page,
                                "select": "id,doi"}).get("results", [])


def doi_str(x):
    if not x:
        return ""
    if isinstance(x, str) and x.lower().startswith("https://doi.org/"):
        return x.split("https://doi.org/")[-1].strip()
    return x.strip()


def slugify(text):
    return re.sub(r'[^A-Za-z0-9._-]+', '_', text)


# ============================================================
# Logging (console + per-set log file)
# ============================================================

_log_file = None

def log(msg):
    ts = time.strftime("%H:%M:%S")
    line = f"[{ts}] {msg}"
    print(line, flush=True)
    if _log_file:
        _log_file.write(line + "\n")
        _log_file.flush()


# ============================================================
# Drawing (same as notebook version)
# ============================================================

def _prepare_draw_graph(GX):
    if GX.number_of_nodes() <= MAX_NODES_TO_DRAW:
        return GX
    keep = set(n for n, d in GX.nodes(data=True) if d.get("type") == "focal")
    depth_sorted = sorted([(n, GX.nodes[n].get("depth", 99)) for n in GX.nodes() if n not in keep],
                          key=lambda x: x[1])
    for n, _ in depth_sorted[:max(0, MAX_NODES_TO_DRAW - len(keep))]:
        keep.add(n)
    return GX.subgraph(keep).copy()


def draw_with_options(GX, title, path, color_by_depth=False):
    GX = _prepare_draw_graph(GX)
    if GX.number_of_nodes() < 1200:
        pos = nx.spring_layout(GX, seed=42, k=0.7)
    else:
        # spectral_layout reliably fails to converge (ArpackNoConvergence) on this data's
        # graph structure, costing ~40-50s per attempt before falling back anyway — skip
        # straight to the fast, reliable fallback instead of paying that cost every time.
        pos = nx.spring_layout(GX, seed=42, k=0.5, iterations=20)

    plt.figure(figsize=(12, 10))
    focal_nodes = [n for n, d in GX.nodes(data=True) if d.get("type") == "focal"]
    other_nodes = [n for n in GX.nodes() if n not in focal_nodes]

    node_colors = None
    legend_handles = []
    if color_by_depth:
        depths = sorted({GX.nodes[n].get("depth") for n in other_nodes if GX.nodes[n].get("depth") is not None})
        cmap = plt.get_cmap("tab10")
        depth_to_color = {d: cmap(i % 10) for i, d in enumerate(depths)}
        node_colors = [depth_to_color.get(GX.nodes[n].get("depth"), "#CCCCCC") for n in other_nodes]
        for d in depths:
            legend_handles.append(Patch(facecolor=depth_to_color[d], edgecolor="none", label=f"Depth {d}"))

    nx.draw_networkx_nodes(GX, pos, nodelist=other_nodes, node_size=28, node_color=node_colors)
    nx.draw_networkx_edges(GX, pos, arrows=True, arrowstyle="-|>", arrowsize=10, width=0.35)
    nx.draw_networkx_nodes(GX, pos, nodelist=focal_nodes, node_size=650, node_shape="s")
    labels = {n: (GX.nodes[n].get("doi") or n.split("/")[-1]) for n in focal_nodes}
    nx.draw_networkx_labels(GX, pos, labels=labels, font_size=8)

    legend_handles.insert(0, Line2D([0], [0], marker="s", linestyle="None", markersize=10, label="Source DOI"))
    legend_handles.append(Line2D([0], [0], marker="o", linestyle="None", markersize=6, label="Alters"))
    plt.legend(handles=legend_handles, loc="lower left", fontsize=8, frameon=False)

    plt.title(title)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(path, dpi=170)
    plt.close()


# ============================================================
# Checkpointing
# ============================================================

def save_progress_csvs(out_dir, set_label, G, depth_by_node):
    depth_counts = collections.Counter(depth_by_node.values())
    pd.DataFrame([{"depth": d, "count": c} for d, c in sorted(depth_counts.items())]).to_csv(
        os.path.join(out_dir, f"{set_label}_nodes_by_depth_IN.csv"), index=False
    )
    # heapq.nlargest is O(n log k) with k=50, vs. the old full sort which was O(n log n) —
    # on an 80,000+ node graph, checkpointing every 25 targets, the old approach meant
    # repeatedly sorting the ENTIRE graph from scratch hundreds of times over a run.
    # This was almost certainly the real cause of the memory/CPU blowup on SetB.
    top50_ids = heapq.nlargest(50, G.in_degree(), key=lambda x: x[1])
    rows = [{"node": n, "in_degree": d, "type": G.nodes[n].get("type", ""), "doi": G.nodes[n].get("doi", "")}
            for n, d in top50_ids]
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, f"{set_label}_top_hubs_IN.csv"), index=False)


# ============================================================
# Main expansion routine (multi-threaded, progress-reporting)
# ============================================================

def save_progress_csvs_streaming(out_dir, set_label, depth_by_node, in_degree):
    """Lightweight checkpoint from in-memory dicts only (no full graph needed)."""
    depth_counts = collections.Counter(depth_by_node.values())
    pd.DataFrame([{"depth": d, "count": c} for d, c in sorted(depth_counts.items())]).to_csv(
        os.path.join(out_dir, f"{set_label}_nodes_by_depth_IN.csv"), index=False
    )
    top50 = heapq.nlargest(50, in_degree.items(), key=lambda x: x[1])
    pd.DataFrame([{"node": n, "in_degree": d} for n, d in top50]).to_csv(
        os.path.join(out_dir, f"{set_label}_top_hubs_IN.csv"), index=False
    )


def finalize_top_hubs(out_dir, set_label, in_degree, node_type_doi):
    """Post-fetch: annotate the top-50 hubs with type/doi by looking them up in the
    small focal/recent-node cache we kept, falling back to blank if not cached."""
    top50 = heapq.nlargest(50, in_degree.items(), key=lambda x: x[1])
    rows = []
    for n, d in top50:
        t, doi = node_type_doi.get(n, ("", ""))
        rows.append({"node": n, "in_degree": d, "type": t, "doi": doi})
    pd.DataFrame(rows).to_csv(os.path.join(out_dir, f"{set_label}_top_hubs_IN.csv"), index=False)


def build_sampled_graph_from_files(out_dir, set_label, max_nodes=1200):
    """Build a small networkx graph for plotting by SAMPLING from the on-disk edge/node
    files, rather than ever materializing the full graph in memory. Used only after the
    fetch phase (and its data writes) are already complete and safe."""
    nodes_path = os.path.join(out_dir, f"{set_label}_nodes_attributes_IN.csv")
    edges_path = os.path.join(out_dir, f"{set_label}_edges_IN.csv")

    nodes_df = pd.read_csv(nodes_path)
    if len(nodes_df) > max_nodes:
        focal = nodes_df[nodes_df["type"] == "focal"]
        others = nodes_df[nodes_df["type"] != "focal"].sort_values("depth")
        keep_n = max(0, max_nodes - len(focal))
        sampled = pd.concat([focal, others.head(keep_n)])
    else:
        sampled = nodes_df
    keep_ids = set(sampled["node"])

    GX = nx.DiGraph()
    for _, row in sampled.iterrows():
        GX.add_node(row["node"], type=row["type"], doi=row.get("doi", ""), depth=row.get("depth", 0))

    # Stream the edges file rather than loading it all — only keep edges within the sample
    with open(edges_path) as f:
        reader = csv.reader(f)
        next(reader, None)  # header
        for src, tgt in reader:
            if src in keep_ids and tgt in keep_ids:
                GX.add_edge(src, tgt)

    return GX


def expand_incoming_multihop(set_label, doi_list, in_depth):
    """Streaming version: nodes/edges are written directly to disk as they're discovered.
    Only lightweight ID -> depth/in-degree tracking lives in RAM, regardless of how large
    the graph grows. This replaces the earlier in-memory-networkx-graph approach, which
    hit a hard memory ceiling around ~85,000 nodes on a 7.7GB machine (two separate
    incidents: a slow swap-thrashing stall, and a crash)."""
    global _log_file
    out_dir = os.path.join(OUT_BASE_DIR_IN, set_label)
    ego_dir = os.path.join(out_dir, f"ego_IN_depth{in_depth}")
    os.makedirs(out_dir, exist_ok=True)
    os.makedirs(ego_dir, exist_ok=True)

    log_path = os.path.join(out_dir, f"run_log_{set_label}.txt")
    _log_file = open(log_path, "a")

    edges_path = os.path.join(out_dir, f"{set_label}_edges_IN.csv")
    nodes_path = os.path.join(out_dir, f"{set_label}_nodes_attributes_IN.csv")
    edges_file = open(edges_path, "w", newline="")
    edges_writer = csv.writer(edges_file)
    edges_writer.writerow(["source", "target"])
    nodes_file = open(nodes_path, "w", newline="")
    nodes_writer = csv.writer(nodes_file)
    nodes_writer.writerow(["node", "type", "doi", "depth"])

    depth_by_node = {}            # id -> depth (lightweight: ints, not full node objects)
    in_degree = collections.Counter()   # id -> in-degree (lightweight: ints)
    node_type_doi = {}            # small cache for focal + a bounded number of others, for labeling
    focal_ids = []
    total_edges_count = 0
    write_lock = Lock()

    log(f"[{set_label}] Resolving {len(doi_list)} seed DOIs...")
    for i, doi in enumerate(doi_list, 1):
        try:
            w = work_from_doi(doi)
        except Exception as e:
            log(f"[{set_label}] Failed seed DOI {doi}: {e}")
            continue
        fid = w.get("id", "")
        if fid in depth_by_node:
            continue  # duplicate seed DOI resolving to the same work
        focal_ids.append(fid)
        fdoi = doi_str(w.get("doi")) or doi
        depth_by_node[fid] = 0
        node_type_doi[fid] = ("focal", fdoi)
        nodes_writer.writerow([fid, "focal", fdoi, 0])
        if i % 10 == 0 or i == len(doi_list):
            log(f"[{set_label}] resolved {i}/{len(doi_list)} seeds")
    nodes_file.flush()

    log(f"[{set_label}] Seed resolution complete: {len(focal_ids)}/{len(doi_list)} resolved.")

    frontier = list(focal_ids)
    visited = set()
    current_depth = 0
    run_start = time.time()
    total_node_count = len(depth_by_node)

    while current_depth < in_depth and total_node_count < NODES_HARD_CAP_PER_SET:
        depth_start = time.time()
        targets = [t for t in frontier if t not in visited]
        visited.update(targets)
        effective_checkpoint_every = max(CHECKPOINT_EVERY, len(targets) // 40 or 1)
        log(f"[{set_label}] === Starting depth {current_depth + 1}/{in_depth}: "
            f"{len(targets)} targets to process (workers={MAX_WORKERS}, "
            f"checkpointing every {effective_checkpoint_every}, streaming to disk) ===")

        next_frontier = []
        calls_done = 0
        calls_lock = Lock()
        depth_call_start = time.time()

        def fetch_one(target):
            try:
                return target, get_citers(target, per_page=MAX_CITERS_PER_NODE), None
            except Exception as e:
                return target, [], e

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {executor.submit(fetch_one, t): t for t in targets}
            for future in as_completed(futures):
                if total_node_count >= NODES_HARD_CAP_PER_SET:
                    break

                target, citers, err = future.result()
                del futures[future]  # release the completed future's retained result (raw
                                      # API response data) — without this, the futures dict
                                      # keeps every response alive in memory for the entire
                                      # depth's duration, defeating the point of streaming
                                      # writes to disk. as_completed() already took its own
                                      # internal snapshot, so mutating this dict is safe.
                if err:
                    log(f"[{set_label}] Citers fetch failed for {target}: {err}")

                with write_lock:
                    for cw in citers:
                        cid = cw.get("id", "")
                        if not cid:
                            continue
                        if cid not in depth_by_node:
                            depth_by_node[cid] = current_depth + 1
                            cdoi = doi_str(cw.get("doi")) or ""
                            ctype = f"in_level{current_depth + 1}"
                            nodes_writer.writerow([cid, ctype, cdoi, current_depth + 1])
                            total_node_count += 1
                            # keep a bounded cache for hub labeling later (only nodes likely
                            # to matter — cheap to keep all of them here too, since this is
                            # just two short strings per node, far lighter than a full graph)
                            node_type_doi[cid] = (ctype, cdoi)
                        edges_writer.writerow([cid, target])
                        total_edges_count += 1
                        in_degree[target] += 1
                        if current_depth + 1 < in_depth:
                            next_frontier.append(cid)
                    node_count = total_node_count
                    edge_count = total_edges_count

                with calls_lock:
                    calls_done += 1
                    n = calls_done

                if n % PROGRESS_EVERY == 0 or n == len(targets):
                    elapsed = time.time() - depth_call_start
                    rate = n / elapsed if elapsed > 0 else 0
                    pct = 100 * n / len(targets) if targets else 100
                    eta_sec = (len(targets) - n) / rate if rate > 0 else float("nan")
                    total_elapsed = time.time() - run_start
                    log(f"[{set_label}] depth {current_depth + 1}: {n}/{len(targets)} ({pct:.1f}%) | "
                        f"{rate:.2f} calls/sec | ETA this depth: {eta_sec / 60:.1f} min | "
                        f"nodes={node_count} edges={edge_count} | "
                        f"total elapsed: {total_elapsed / 60:.1f} min")

                if n % effective_checkpoint_every == 0:
                    with write_lock:
                        edges_file.flush()
                        nodes_file.flush()
                        save_progress_csvs_streaming(out_dir, set_label, depth_by_node, in_degree)
                    log(f"[{set_label}] checkpoint saved ({n} targets this depth) — "
                        f"edges/nodes files flushed to disk")

        depth_elapsed = time.time() - depth_start
        log(f"[{set_label}] === Depth {current_depth + 1} complete in {depth_elapsed / 60:.1f} min | "
            f"nodes={total_node_count} edges={total_edges_count} ===")

        frontier = list(set(next_frontier))
        current_depth += 1

        if total_node_count >= NODES_HARD_CAP_PER_SET:
            log(f"[{set_label}] Reached node hard cap ({NODES_HARD_CAP_PER_SET}); stopping IN expansion.")
            break

    edges_file.flush()
    edges_file.close()
    nodes_file.flush()
    nodes_file.close()

    total_elapsed = time.time() - run_start
    log(f"[{set_label}] All depths complete in {total_elapsed / 60:.1f} min total. "
        f"Final: nodes={total_node_count} edges={total_edges_count}")

    # ============================================================
    # Summary + hub CSVs — cheap, from lightweight dicts. Edges/node-attributes files
    # are ALREADY complete on disk (written incrementally throughout the run above) —
    # nothing left to do for data safety at this point.
    # ============================================================
    save_progress_csvs_streaming(out_dir, set_label, depth_by_node, in_degree)
    finalize_top_hubs(out_dir, set_label, in_degree, node_type_doi)
    pd.DataFrame([{
        "set": set_label,
        "direction": "IN",
        "depth": in_depth,
        "focal_count": len(focal_ids),
        "total_nodes": total_node_count,
        "total_edges": total_edges_count,
        "hard_cap_reached": int(total_node_count >= NODES_HARD_CAP_PER_SET)
    }]).to_csv(os.path.join(out_dir, f"{set_label}_summary_IN.csv"), index=False)
    log(f"[{set_label}] Data files finalized (summary, edges, node attributes, top hubs) — "
        f"all safe on disk regardless of plotting outcome below.")

    # ============================================================
    # Plotting — built from a SAMPLE read off disk, never the full in-memory graph.
    # A single combined plot only (per-focal "ego" plots were found to be near-duplicates
    # of the combined plot at this scale and mostly wasted time — skipped for large graphs).
    # ============================================================
    try:
        log(f"[{set_label}] Building sampled graph for plotting (max 1200 nodes)...")
        G_sample = build_sampled_graph_from_files(out_dir, set_label, max_nodes=MAX_NODES_TO_DRAW)
        combined_png = os.path.join(out_dir, f"{set_label}_combined_IN_depth{in_depth}.png")
        draw_with_options(G_sample, f"{set_label} — IN up to depth {in_depth} (sampled)",
                           combined_png, color_by_depth=COLOR_BY_DEPTH)
        log(f"[{set_label}] Combined plot saved (from sample of {G_sample.number_of_nodes()} nodes).")
    except Exception as e:
        log(f"  [{set_label}] combined plot FAILED ({type(e).__name__}: {e}) — skipping, data unaffected")

    if total_node_count <= 3000:
        log(f"[{set_label}] Drawing {len(focal_ids)} per-source ego plots (graph small enough)...")
        ego_failures = 0
        for i, fid in enumerate(focal_ids, 1):
            keep = {fid} | {n for n, d in depth_by_node.items() if d is not None and 1 <= d <= in_depth}
            # only feasible at this size; still sample from files for consistency
            try:
                G_full_small = build_sampled_graph_from_files(out_dir, set_label, max_nodes=total_node_count)
                G_ego = G_full_small.subgraph(keep & set(G_full_small.nodes())).copy()
                t, doi = node_type_doi.get(fid, ("", ""))
                ego_png = os.path.join(ego_dir, f"ego_IN_depth{in_depth}_{slugify(doi or fid.split('/')[-1])}.png")
                draw_with_options(G_ego, f"IN ego (depth {in_depth}) — {doi}", ego_png, color_by_depth=COLOR_BY_DEPTH)
            except Exception as e:
                ego_failures += 1
                log(f"  [{set_label}] ego plot {i}/{len(focal_ids)} FAILED ({type(e).__name__}: {e}) — skipping")
        if ego_failures:
            log(f"[{set_label}] {ego_failures}/{len(focal_ids)} ego plots failed and were skipped")
    else:
        log(f"[{set_label}] Skipping per-source ego plots — graph too large ({total_node_count} nodes) "
            f"for 50 individual renders to be worthwhile; combined sampled plot above is representative.")

    log(f"[{set_label}] IN: Saved combined map + CSVs at {out_dir}")

    _log_file.close()
    _log_file = None

    return focal_ids, depth_by_node


# ============================================================
# CLI entry point
# ============================================================

def run_set(label, dois, in_depth):
    log(f"\n{'=' * 60}\nStarting {label} (n={len(dois)} seeds, depth={in_depth})\n{'=' * 60}")
    result = expand_incoming_multihop(label, dois, in_depth)
    zip_name = f"{label}_results"
    shutil.make_archive(zip_name, "zip", os.path.join(OUT_BASE_DIR_IN, label))
    log(f"[{label}] Zipped results to {zip_name}.zip")
    return result


def main():
    global MAX_WORKERS

    parser = argparse.ArgumentParser(description="Run the IN-only multi-hop citation map builder locally.")
    parser.add_argument("--set", choices=["A", "B", "both"], default="both",
                        help="Which set to run (default: both, sequentially)")
    parser.add_argument("--depth", type=int, default=IN_DEPTH_DEFAULT,
                        help=f"IN_DEPTH override (default: {IN_DEPTH_DEFAULT})")
    parser.add_argument("--workers", type=int, default=MAX_WORKERS,
                        help=f"Concurrent request workers (default: {MAX_WORKERS})")
    args = parser.parse_args()

    MAX_WORKERS = args.workers

    os.makedirs(OUT_BASE_DIR_IN, exist_ok=True)

    if CONTACT_EMAIL == "your_email@nccu.edu":
        print("!! WARNING: CONTACT_EMAIL is still the placeholder. Edit the script and set your real "
              "email before running a long job — OpenAlex's polite pool gives much better reliability. !!")
        time.sleep(3)

    if args.set in ("A", "both"):
        run_set(LABEL_A, DOIS_SET_A, args.depth)
    if args.set in ("B", "both"):
        run_set(LABEL_B, DOIS_SET_B, args.depth)

    print("\n✅ All requested runs complete.")


if __name__ == "__main__":
    main()
