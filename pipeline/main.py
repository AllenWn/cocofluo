#!/usr/bin/env python3
"""
Cocofluo Pipeline - Main Entry Point

Batch pipeline:
1) Read multi-channel confocal TIFF
2) Extract structural (DAPI) channel -> optional contrast enhancement
3) Cell segmentation via Cellpose (instance mask)
4) For each marker channel: percentile threshold -> binary marker mask
5) For each marker: compute "express" cells by ratio(marker_area_in_cell / cell_area) >= threshold
6) Colocalization: intersection counts for configured marker groups
7) Write per-image report.csv
8) Sampling QC images every N files:
   - cellseg: (left raw DAPI, right colored segmentation overlay on raw)
   - markerextract: (left raw marker, right colored binary mask)
   - express: (left raw DAPI, right highlight express cells on raw)
"""

import os
import sys
import time
import glob
import csv
import colorsys
import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cellpose import models, io

import config as cfg


def ch_to_idx(ch_1based: int) -> int:
    """1-based channel -> 0-based index."""
    return int(ch_1based) - 1


def get_channel_name(ch_1based: int) -> str:
    """Human-readable channel name from config; fallback to chK."""
    name = cfg.CHANNEL_NAMES.get(int(ch_1based), "")
    return name if name else f"ch{int(ch_1based)}"


def normalize_channel(img2d, p_low=1, p_high=99):
    """Robust percentile normalization -> float32 in [0,1]."""
    arr = img2d.astype(np.float32)
    p1, p99 = np.percentile(arr, (p_low, p_high))
    out = np.clip((arr - p1) / (p99 - p1 + 1e-8), 0.0, 1.0)
    return out.astype(np.float32)


def safe_minmax01(img2d):
    """Fallback normalization -> float32 in [0,1]."""
    x = img2d.astype(np.float32)
    mn, mx = float(np.min(x)), float(np.max(x))
    if mx <= mn + 1e-8:
        return np.zeros_like(x, dtype=np.float32)
    return ((x - mn) / (mx - mn)).astype(np.float32)


def extract_marker_binary(marker2d, percentile):
    """
    Percentile thresholding on marker channel.
    Returns (binary uint8 mask 0/1, threshold float).
    """
    arr = marker2d.astype(np.float32)
    thr = float(np.percentile(arr, float(percentile)))
    binmask = (arr > thr).astype(np.uint8)
    return binmask, thr


def compute_express_ids(instance_mask, marker_bin, ratio_threshold):
    """
    For each cell id, compute marker-positive area ratio:
        ratio = sum(marker_bin in cell) / cell_area
    Express if ratio >= ratio_threshold

    Returns:
        express_ids (set[int])
        ratios (dict[int,float])  # for debugging
    """
    masks = instance_mask
    marker = marker_bin

    # bincount index = cell_id; cell_id=0 is background
    cell_area = np.bincount(masks.ravel())
    marker_in_cell = np.bincount(
        masks.ravel(), weights=marker.astype(np.float64).ravel()
    )

    cell_ids = np.unique(masks)
    cell_ids = cell_ids[cell_ids != 0]

    express = set()
    ratios = {}
    for cid in cell_ids:
        cid = int(cid)
        a = cell_area[cid] if cid < len(cell_area) else 0
        if a <= 0:
            continue
        r = float(marker_in_cell[cid] / a)
        ratios[cid] = r
        if r >= float(ratio_threshold):
            express.add(cid)
    return express, ratios


def split_channels(img):
    """
    Return (n_channels, get_channel_func)

    Supports common layouts:
    - (C, H, W)
    - (H, W, C)
    Heuristic:
      - if last dim is 3/4 and first dim is large -> (H,W,C)
      - else if first dim is small (<=32) and smaller than H/W -> (C,H,W)
      - else -> (H,W,C)
    """
    if img.ndim == 2:
        raise ValueError("Expected multi-channel image; got 2D.")
    if img.ndim != 3:
        raise ValueError(
            f"Unsupported image ndim={img.ndim}, shape={img.shape}. Expected 3D."
        )

    # prefer (H,W,C) if last dim looks like RGB(A) and first dim not tiny
    if img.shape[-1] in (3, 4) and img.shape[0] > 4:

        def get_ch(idx0):
            return img[:, :, idx0]

        return img.shape[2], get_ch

    # prefer (C,H,W) if first dim looks like channels
    if (
        (img.shape[0] <= 32)
        and (img.shape[0] < img.shape[1])
        and (img.shape[0] < img.shape[2])
    ):

        def get_ch(idx0):
            return img[idx0, :, :]

        return img.shape[0], get_ch

    # fallback: treat as (H,W,C)
    def get_ch(idx0):
        return img[:, :, idx0]

    return img.shape[2], get_ch


def hsv_to_rgb_array(h, s, v):
    h = np.asarray(h)
    s = np.asarray(s)
    v = np.asarray(v)
    out = np.zeros(h.shape + (3,), dtype=np.float32)
    flat_h, flat_s, flat_v = h.ravel(), s.ravel(), v.ravel()
    out_flat = out.reshape(-1, 3)
    for i in range(flat_h.size):
        out_flat[i] = colorsys.hsv_to_rgb(
            float(flat_h[i]), float(flat_s[i]), float(flat_v[i])
        )
    return (out * 255.0).clip(0, 255).astype(np.uint8)


def make_bright_mask_rgb(masks, seed=123, s=1.0, v=1.0):
    """
    Instance mask -> bright RGB (no white).
    Background = black.
    """
    n = int(masks.max())
    rng = np.random.RandomState(seed)
    hues = rng.rand(n + 1).astype(np.float32)
    hues[0] = 0.0  # bg
    hsv_s = np.full_like(hues, fill_value=float(s), dtype=np.float32)
    hsv_v = np.full_like(hues, fill_value=float(v), dtype=np.float32)
    lut = hsv_to_rgb_array(hues, hsv_s, hsv_v)  # (n+1,3)
    rgb = lut[masks]  # (H,W,3)
    rgb[masks == 0] = 0
    return rgb


def make_overlay(gray01, mask_rgb, alpha=0.55):
    """Overlay colored mask on grayscale base (gray01 in [0,1])."""
    gray_u8 = (gray01 * 255.0).clip(0, 255).astype(np.uint8)
    gray_rgb = np.stack([gray_u8, gray_u8, gray_u8], axis=-1).astype(np.float32)

    mask_f = mask_rgb.astype(np.float32)
    m2d = (mask_rgb.sum(axis=-1) > 0).astype(np.float32)  # (H,W)
    m3d = m2d[..., None]  # (H,W,1)

    out = gray_rgb * (1.0 - float(alpha) * m3d) + mask_f * (float(alpha) * m3d)
    return out.clip(0, 255).astype(np.uint8)


def to_u8(img2d):
    """Linear stretch to uint8 for visualization."""
    x = img2d.astype(np.float32)
    mn, mx = float(np.min(x)), float(np.max(x))
    if mx <= mn + 1e-8:
        return np.zeros_like(x, dtype=np.uint8)
    y = (x - mn) / (mx - mn)
    return (y * 255.0).clip(0, 255).astype(np.uint8)


def make_color_mask(bin_mask, color=(255, 0, 255)):
    """Binary mask -> solid-color RGB image."""
    rgb = np.zeros((bin_mask.shape[0], bin_mask.shape[1], 3), dtype=np.uint8)
    rgb[bin_mask > 0] = np.array(color, dtype=np.uint8)
    return rgb


def colorize_express_cells(instance_mask, express_ids, seed=123):
    """Only express cells colored; others black."""
    H, W = instance_mask.shape
    out = np.zeros((H, W, 3), dtype=np.uint8)
    if not express_ids:
        return out

    rng = np.random.default_rng(seed)
    colors = {}
    for cid in express_ids:
        c = rng.integers(40, 256, size=3, dtype=np.uint8)
        # prevent near-white: force one channel low
        k = int(rng.integers(0, 3))
        c[k] = rng.integers(0, 80, dtype=np.uint8)
        colors[int(cid)] = c

    for cid in express_ids:
        out[instance_mask == int(cid)] = colors[int(cid)]
    return out


def save_cell_seg_sampling(cell_raw, masks, base_name, out_dir):
    """
    One PNG with 2 panels:
      Left: raw DAPI (grayscale)
      Right: segmentation result (colored instances overlaid on raw)
    """
    raw_u8 = to_u8(cell_raw)
    raw01 = raw_u8.astype(np.float32) / 255.0

    mask_rgb = make_bright_mask_rgb(masks, seed=123, s=1.0, v=1.0)
    overlay = make_overlay(raw01, mask_rgb, alpha=0.55)

    plt.figure(figsize=(18, 8))
    plt.subplot(1, 2, 1)
    plt.imshow(raw_u8, cmap="gray")
    plt.axis("off")

    plt.subplot(1, 2, 2)
    plt.imshow(overlay)
    plt.axis("off")

    p = os.path.join(out_dir, f"{base_name}_cellseg_sampling.png")
    plt.savefig(p, dpi=150, bbox_inches="tight", pad_inches=0)
    plt.close()
    return p


def save_marker_extract_sampling(
    marker_raw, marker_bin, marker_ch, thr, base_name, out_dir
):
    """
    One PNG with 2 panels:
      Left: raw marker (grayscale)
      Right: extracted binary mask colored (NO overlay)
    """
    raw_u8 = to_u8(marker_raw)
    mask_rgb = make_color_mask(marker_bin, color=(255, 0, 255))  # bright magenta

    plt.figure(figsize=(18, 7))
    plt.subplot(1, 2, 1)
    plt.imshow(raw_u8, cmap="gray")
    plt.axis("off")

    plt.subplot(1, 2, 2)
    plt.imshow(mask_rgb)
    plt.axis("off")

    ch_name = get_channel_name(marker_ch)
    p = os.path.join(out_dir, f"{base_name}_markerextract_{ch_name}_sampling.png")
    plt.savefig(p, dpi=150, bbox_inches="tight", pad_inches=0)
    plt.close()
    return p


def save_express_sampling(cell_raw, masks, express_ids, marker_ch, base_name, out_dir):
    """
    One PNG with 2 panels:
      Left: raw DAPI (grayscale)
      Right: only express cells colored on top of DAPI (others unchanged)
    """
    raw_u8 = to_u8(cell_raw)
    base_rgb = np.stack([raw_u8, raw_u8, raw_u8], axis=-1)

    express_rgb = colorize_express_cells(masks, express_ids, seed=123)
    region = express_rgb.sum(axis=-1) > 0

    result_rgb = base_rgb.copy()
    result_rgb[region] = express_rgb[region]

    plt.figure(figsize=(18, 8))
    plt.subplot(1, 2, 1)
    plt.imshow(raw_u8, cmap="gray")
    plt.axis("off")

    plt.subplot(1, 2, 2)
    plt.imshow(result_rgb)
    plt.axis("off")

    ch_name = get_channel_name(marker_ch)
    p = os.path.join(out_dir, f"{base_name}_express_{ch_name}_sampling.png")
    plt.savefig(p, dpi=150, bbox_inches="tight", pad_inches=0)
    plt.close()
    return p


def build_report_header(marker_channels_1b, coloc_groups):
    headers = ["file", "num_cells"]
    for ch in marker_channels_1b:
        headers.append(f"express_{get_channel_name(ch)}")
    for group in coloc_groups:
        names = [get_channel_name(ch) for ch in group]
        headers.append("coloc_" + "_".join(names))
    return headers


def build_report_row(result, marker_channels_1b, coloc_groups):
    row = [result["file"], result["num_cells"]]
    for ch in marker_channels_1b:
        row.append(int(result["express_counts"].get(ch, 0)))
    for group in coloc_groups:
        row.append(int(result["coloc_counts"].get(tuple(group), 0)))
    return row


def choose_sampling_marker(marker_channels_1b):
    """
    Pick a marker channel for sampling:
    prefer channels with non-empty names, otherwise first available.
    """
    named = [ch for ch in marker_channels_1b if cfg.CHANNEL_NAMES.get(int(ch), "")]
    return (
        named[0] if named else (marker_channels_1b[0] if marker_channels_1b else None)
    )


def sanitize_coloc_groups(groups, dapi_ch, valid_channels_set):
    """
    Drop invalid coloc groups:
      - contains DAPI
      - contains channel not in valid set
      - duplicated channel in group
    Return list[tuple[int,...]]
    """
    out = []
    for g in groups:
        g = tuple(int(x) for x in g)
        if int(dapi_ch) in g:
            print(f"[WARN] Skip coloc group {g}: contains DAPI channel {dapi_ch}")
            continue
        if any((ch not in valid_channels_set) for ch in g):
            print(f"[WARN] Skip coloc group {g}: channel out of range")
            continue
        if len(set(g)) != len(g):
            print(f"[WARN] Skip coloc group {g}: duplicated channel")
            continue
        out.append(g)
    return out


def process_single_image(image_path, model, marker_channels_1b, coloc_groups):
    """
    Process one TIFF. Returns dict:
      - file, num_cells
      - express_counts: {ch:count}
      - coloc_counts: {tuple(chs):count}
    And payload for sampling.
    """
    base_name = os.path.splitext(os.path.basename(image_path))[0]

    img = io.imread(image_path)
    n_ch_in_file, get_ch = split_channels(img)

    # DAPI channel
    dapi_idx0 = ch_to_idx(cfg.DAPI_CH)
    if dapi_idx0 < 0 or dapi_idx0 >= n_ch_in_file:
        raise ValueError(
            f"DAPI_CH={cfg.DAPI_CH} out of range for file channels={n_ch_in_file}"
        )

    cell_raw = get_ch(dapi_idx0).astype(np.float32)

    if bool(cfg.DAPI_ENHANCEMENT):
        cell_norm = normalize_channel(cell_raw, 1, 99)
    else:
        cell_norm = safe_minmax01(cell_raw)

    # Cellpose v4 returns 3 items: masks, flows, styles
    masks, flows, styles = model.eval(
        cell_norm,
        diameter=None,
        flow_threshold=float(cfg.FLOW_THRESHOLD),
        cellprob_threshold=float(cfg.CELLPROB_THRESHOLD),
    )

    num_cells = int(masks.max())

    # Marker extraction + express
    marker_bins = {}
    marker_thrs = {}
    marker_raws = {}
    express_sets = {}

    for ch_1b in marker_channels_1b:
        idx0 = ch_to_idx(ch_1b)
        if idx0 < 0 or idx0 >= n_ch_in_file:
            continue

        marker_raw = get_ch(idx0).astype(np.float32)

        binmask, thr = extract_marker_binary(
            marker_raw, float(cfg.MARKER_EXTRACT_PERCENTILE)
        )
        marker_bins[int(ch_1b)] = binmask
        marker_thrs[int(ch_1b)] = thr
        marker_raws[int(ch_1b)] = marker_raw

        express_ids, _ = compute_express_ids(
            masks, binmask, float(cfg.MARKER_EXPRESS_THRESHOLD)
        )
        express_sets[int(ch_1b)] = express_ids

    express_counts = {
        int(ch): len(express_sets.get(int(ch), set())) for ch in marker_channels_1b
    }

    # Colocalization: intersection
    coloc_counts = {}
    for group in coloc_groups:
        group = tuple(int(x) for x in group)
        if not all((ch in express_sets) for ch in group):
            coloc_counts[group] = 0
            continue
        common = set(express_sets[group[0]])
        for ch in group[1:]:
            common &= express_sets[ch]
        coloc_counts[group] = len(common)

    return {
        "file": base_name,
        "num_cells": num_cells,
        "express_counts": express_counts,
        "coloc_counts": coloc_counts,
        # sampling payload
        "_cell_raw": cell_raw,
        "_cell_norm": cell_norm,
        "_masks": masks,
        "_marker_bins": marker_bins,
        "_marker_thrs": marker_thrs,
        "_express_sets": express_sets,
        "_marker_raws": marker_raws,
        "_n_ch_in_file": n_ch_in_file,
    }


def main():
    print("=" * 60)
    print("Cocofluo Pipeline")
    print("=" * 60)

    os.makedirs(cfg.OUTPUT_DIR, exist_ok=True)
    os.makedirs(cfg.SAMPLE_DIR, exist_ok=True)

    pattern = os.path.join(cfg.INPUT_DIR, "*.tif")
    image_files = sorted(glob.glob(pattern))
    if not image_files:
        print(f"[ERROR] No .tif found in: {cfg.INPUT_DIR}")
        sys.exit(1)

    print(f"[INFO] Found {len(image_files)} images in {cfg.INPUT_DIR}")

    # --- Determine effective channels from the FIRST file (avoid cfg.N_CHANNELS mismatch) ---
    first_img = io.imread(image_files[0])
    n_ch_first, _ = split_channels(first_img)

    if int(cfg.N_CHANNELS) != int(n_ch_first):
        print(
            f"[WARN] cfg.N_CHANNELS={cfg.N_CHANNELS} but first file has {n_ch_first} channels. "
            f"Using file channel count for marker list/report columns."
        )

    dapi_ch = int(cfg.DAPI_CH)
    if dapi_ch < 1 or dapi_ch > n_ch_first:
        print(
            f"[ERROR] cfg.DAPI_CH={cfg.DAPI_CH} out of range for first file channels={n_ch_first}"
        )
        sys.exit(1)

    # Marker channels: based on actual file channel count (1-based), excluding DAPI
    marker_channels = [ch for ch in range(1, n_ch_first + 1) if ch != dapi_ch]

    # Sanitize coloc groups (based on valid channels in file)
    valid_channels_set = set(range(1, n_ch_first + 1))
    coloc_groups = sanitize_coloc_groups(
        cfg.COLOCALIZATION_GROUPS, dapi_ch, valid_channels_set
    )

    print(f"[INFO] DAPI channel: ch{dapi_ch} ({get_channel_name(dapi_ch)})")
    print(
        f"[INFO] Marker channels: {[f'ch{ch}({get_channel_name(ch)})' for ch in marker_channels]}"
    )
    print(f"[INFO] Colocalization groups (effective): {coloc_groups}")

    # Load model
    print(f"[INFO] Loading Cellpose model (GPU={bool(cfg.USE_GPU)}) ...")
    model = models.CellposeModel(gpu=bool(cfg.USE_GPU))
    print("[INFO] Model loaded.")

    # Report setup
    headers = build_report_header(marker_channels, coloc_groups)
    rows = []

    SAMPLE_EVERY = 10
    sampling_marker = choose_sampling_marker(marker_channels)

    total = len(image_files)
    ok = 0

    for i, path in enumerate(image_files):
        base_name = os.path.splitext(os.path.basename(path))[0]
        print(f"\n[{i+1}/{total}] Processing: {base_name}")
        t0 = time.time()

        try:
            result = process_single_image(path, model, marker_channels, coloc_groups)
        except Exception as e:
            print(f"  [ERROR] {e}")
            continue

        dt = time.time() - t0
        ok += 1
        print(f"  Cells: {result['num_cells']} | Time: {dt:.2f}s")

        # console stats
        for ch in marker_channels:
            print(
                f"  Express {get_channel_name(ch)}: {int(result['express_counts'].get(ch, 0))}"
            )
        for group in coloc_groups:
            names = [get_channel_name(int(x)) for x in group]
            print(
                f"  Coloc {'+'.join(names)}: {int(result['coloc_counts'].get(tuple(group), 0))}"
            )

        rows.append(build_report_row(result, marker_channels, coloc_groups))

        # Sampling QC: every 10 images (1st, 11th, 21st, ...)
        if (i % SAMPLE_EVERY) == 0:
            print("  [SAMPLING] Saving QC images ...")

            if bool(cfg.CELL_SEG_SAMPLING):
                p = save_cell_seg_sampling(
                    result["_cell_raw"], result["_masks"], base_name, cfg.SAMPLE_DIR
                )
                print(f"    -> {p}")

            if bool(cfg.MARKER_EXTRACT_SAMPLING) and sampling_marker is not None:
                ch = int(sampling_marker)
                if ch in result["_marker_bins"]:
                    p = save_marker_extract_sampling(
                        result["_marker_raws"][ch],
                        result["_marker_bins"][ch],
                        ch,
                        float(result["_marker_thrs"][ch]),
                        base_name,
                        cfg.SAMPLE_DIR,
                    )
                    print(f"    -> {p}")

            if bool(cfg.MARKER_EXPRESS_SAMPLING) and sampling_marker is not None:
                ch = int(sampling_marker)
                if ch in result["_express_sets"]:
                    p = save_express_sampling(
                        result["_cell_raw"],
                        result["_masks"],
                        result["_express_sets"][ch],
                        ch,
                        base_name,
                        cfg.SAMPLE_DIR,
                    )
                    print(f"    -> {p}")

    # Write report
    report_path = os.path.join(cfg.OUTPUT_DIR, "report.csv")
    print(f"\n[INFO] Writing report: {report_path}")
    with open(report_path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(headers)
        for r in rows:
            w.writerow(r)

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"Images found: {len(image_files)}")
    print(f"Images processed OK: {ok}")
    print(f"Report rows written: {len(rows)}")

    if rows:
        total_cells = sum(int(r[1]) for r in rows)
        print(f"Total cells detected: {total_cells}")

        # marker totals
        for j, ch in enumerate(marker_channels):
            col_idx = 2 + j
            total_express = sum(int(r[col_idx]) for r in rows)
            print(f"Total express {get_channel_name(ch)}: {total_express}")

        # coloc totals
        coloc_start = 2 + len(marker_channels)
        for k, group in enumerate(coloc_groups):
            col_idx = coloc_start + k
            names = [get_channel_name(int(x)) for x in group]
            total_coloc = sum(int(r[col_idx]) for r in rows)
            print(f"Total coloc {'+'.join(names)}: {total_coloc}")

    print("\n[DONE]")


if __name__ == "__main__":
    main()
