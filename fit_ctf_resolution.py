#!/usr/bin/env python3
import argparse
import math
import re
from dataclasses import dataclass

import numpy as np
import mrcfile


@dataclass
class Ctffind3Params:
    df1_A: float
    df2_A: float
    angast_deg: float


def parse_ctffind3_log(log_path: str) -> Ctffind3Params:
    pat = re.compile(
        r"^\s*([0-9.+-Ee]+)\s+([0-9.+-Ee]+)\s+([0-9.+-Ee]+)\s+([0-9.+-Ee]+)\s+Final Values",
        re.MULTILINE
    )
    txt = open(log_path, "r", errors="ignore").read()
    m = pat.search(txt)
    if not m:
        raise ValueError(f"Could not find 'Final Values' line in {log_path}")
    return Ctffind3Params(df1_A=float(m.group(1)), df2_A=float(m.group(2)), angast_deg=float(m.group(3)))


def electron_wavelength_A(kv: float) -> float:
    V = kv * 1000.0
    return 12.2639 / math.sqrt(V * (1.0 + 0.97845e-6 * V))


def make_freq_grids(n: int, pixel_A: float):
    cy = (n - 1) / 2.0
    cx = (n - 1) / 2.0
    yy, xx = np.indices((n, n), dtype=np.float32)
    dy = yy - cy
    dx = xx - cx
    theta = np.arctan2(dy, dx)
    kx = dx / (n * pixel_A)  # cycles/Å
    ky = dy / (n * pixel_A)
    k2 = kx * kx + ky * ky
    k = np.sqrt(k2)
    return k, theta, k2


def ctf_astigmatic(k2, theta, df1_A, df2_A, angast_rad, cs_A, lam_A, amp_contrast):
    df_avg = 0.5 * (df1_A + df2_A)
    df_diff = 0.5 * (df1_A - df2_A)
    df_theta = df_avg + df_diff * np.cos(2.0 * (theta - angast_rad))

    chi = np.pi * (lam_A * df_theta * k2 + 0.5 * cs_A * (lam_A ** 3) * (k2 ** 2))

    a = float(amp_contrast)
    a = max(min(a, 1.0), -1.0)
    ctf = -(math.sqrt(max(0.0, 1.0 - a * a)) * np.sin(chi) + a * np.cos(chi))
    return ctf


def annulus_corr(obs2d, pred2d, k, kmin, kmax, nbins, bg_smooth_bins=9):
    k_flat = k.ravel()
    obs_flat = obs2d.ravel()

    edges = np.linspace(kmin, kmax, nbins + 1)
    bin_id = np.digitize(k_flat, edges) - 1
    valid = (bin_id >= 0) & (bin_id < nbins)

    # radial mean background
    rad_mean = np.full(nbins, np.nan, dtype=np.float64)
    for i in range(nbins):
        sel = valid & (bin_id == i)
        if np.sum(sel) > 50:
            rad_mean[i] = float(np.mean(obs_flat[sel]))

    # smooth background
    bg = rad_mean.copy()
    half = max(1, bg_smooth_bins // 2)
    for i in range(nbins):
        lo = max(0, i - half)
        hi = min(nbins, i + half + 1)
        w = bg[lo:hi]
        w = w[np.isfinite(w)]
        bg[i] = float(np.mean(w)) if w.size else np.nan

    bg_per_pix = np.full_like(obs_flat, np.nan, dtype=np.float64)
    for i in range(nbins):
        sel = valid & (bin_id == i)
        if np.isfinite(bg[i]) and np.any(sel):
            bg_per_pix[sel] = bg[i]

    obs_w = obs_flat - bg_per_pix
    pred_flat = (pred2d * pred2d).ravel()  # compare to CTF^2

    centers = 0.5 * (edges[:-1] + edges[1:])
    corr = np.full(nbins, np.nan, dtype=np.float64)

    for i in range(nbins):
        sel = valid & (bin_id == i)
        if np.sum(sel) < 200:
            continue
        o = obs_w[sel]
        p = pred_flat[sel]
        o = o - np.mean(o)
        p = p - np.mean(p)
        so = np.std(o)
        sp = np.std(p)
        if so < 1e-12 or sp < 1e-12:
            continue
        corr[i] = float(np.mean((o / so) * (p / sp)))

    return centers, corr


def fit_resolution_from_corr(k_centers, corr, thr=0.3, consec_fail=4):
    last_good = None
    fails = 0
    for kc, cc in zip(k_centers, corr):
        if not np.isfinite(cc):
            continue
        if cc >= thr:
            last_good = kc
            fails = 0
        else:
            fails += 1
            if fails >= consec_fail:
                break
    if last_good is None or last_good <= 0:
        return float("nan")
    return 1.0 / last_good


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--power_mrc", required=True, help="CTFFIND3 diagnostic spectrum MRC (e.g., *_ctf.mrc or *_power.mrc)")
    ap.add_argument("--log", required=True, help="CTFFIND3 log with 'Final Values'")
    ap.add_argument("--pixel_A", type=float, required=True)
    ap.add_argument("--kv", type=float, required=True)
    ap.add_argument("--cs_mm", type=float, required=True)
    ap.add_argument("--amp_contrast", type=float, required=True)
    ap.add_argument("--resmin_A", type=float, default=50.0, help="Lowest resolution used for scoring (Å)")
    ap.add_argument("--nbins", type=int, default=120)
    ap.add_argument("--thr", type=float, default=0.3)
    ap.add_argument("--consec_fail", type=int, default=4)
    ap.add_argument("--out_csv", default=None)
    args = ap.parse_args()

    p = parse_ctffind3_log(args.log)

    with mrcfile.open(args.power_mrc, permissive=True) as m:
        img = np.array(m.data, dtype=np.float32)
    if img.ndim == 3:
        img = img[0]
    if img.shape[0] != img.shape[1]:
        raise ValueError("power_mrc must be square.")
    n = img.shape[0]

    k, theta, k2 = make_freq_grids(n, args.pixel_A)
    nyquist = 0.5 / args.pixel_A
    kmin = 1.0 / args.resmin_A
    kmax = 0.95 * nyquist

    lam = electron_wavelength_A(args.kv)
    cs_A = args.cs_mm * 1e7  # mm -> Å
    ang = np.deg2rad(p.angast_deg)

    ctf = ctf_astigmatic(k2, theta, p.df1_A, p.df2_A, ang, cs_A, lam, args.amp_contrast)

    k_centers, corr = annulus_corr(img, ctf, k, kmin, kmax, args.nbins)
    fit_res = fit_resolution_from_corr(k_centers, corr, args.thr, args.consec_fail)

    print(f"DF1={p.df1_A:.2f} Å  DF2={p.df2_A:.2f} Å  ANGAST={p.angast_deg:.2f}°")
    print(f"CTF-fit resolution (corr>={args.thr}): {fit_res:.2f} Å")

    if args.out_csv:
        import csv
        with open(args.out_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["k_1_per_A", "res_A", "corr"])
            for kc, cc in zip(k_centers, corr):
                w.writerow([kc, (1.0 / kc) if kc > 0 else np.nan, cc])
        print(f"Wrote: {args.out_csv}")


if __name__ == "__main__":
    main()
