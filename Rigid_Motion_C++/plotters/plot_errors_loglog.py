#!/usr/bin/env python3
import argparse, os, csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42

def read_wide_csv(path):
    with open(path, 'r', newline='') as f:
        rdr = csv.DictReader(f)
        hs, v1, v2, v3, v4, v5 = [], [], [], [], [], []
        for row in rdr:
            hs.append(float(row['h']))
            v1.append(float(row['vanilla']))
            v2.append(float(row['vanilla_feedback']))
            v3.append(float(row['feedback']))
            v4.append(float(row['adaptive']))
            v5.append(float(row['strang']))
    return (np.array(hs),
            np.array(v1), np.array(v2), np.array(v3), np.array(v4), np.array(v5))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv', default='metrics/errors_vs_h.csv',
                    help='CSV path produced by main.cpp')
    ap.add_argument('--out', default=None,
                    help='If set, save to this path (extension decides format, e.g., .png/.pdf)')
    ap.add_argument('--title', default='Error vs Step Size',
                    help='Figure title')
    args = ap.parse_args()

    h, e_van, e_vanfb, e_fb, e_ad, e_str = read_wide_csv(args.csv)

    # Filter zeros or negatives just in case
    mask = (h > 0) & (e_van > 0) & (e_vanfb > 0) & (e_fb > 0) & (e_ad > 0) & (e_str > 0)
    h = h[mask]
    e_van   = e_van[mask]
    e_vanfb = e_vanfb[mask]
    e_fb    = e_fb[mask]
    e_ad    = e_ad[mask]
    e_str   = e_str[mask]

    fig, ax = plt.subplots(figsize=(12,8), constrained_layout=True)
    plt.loglog(h, e_van,   marker='o', markersize=12, linewidth=2.4, label='Euler\'s method')
    plt.loglog(h, e_vanfb, marker='s', markersize=12, linewidth=2.4, label='Feedback (unity gain)')
    plt.loglog(h, e_fb,    marker='^', markersize=12, linewidth=2.4, label='Feedback $(1/hL)$')
    plt.loglog(h, e_ad,    marker='D', markersize=12, linewidth=2.4, label='Adaptive Feedback')
    plt.loglog(h, e_str,   marker='v', markersize=12, linewidth=2.4, label='Splitting method')

    plt.xlabel('$h$', fontsize=26)
    plt.ylabel('max $V(x_k)$', fontsize=26)
    # plt.title(args.title)
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.grid(True, which='both', linestyle='--', alpha=0.8)
    plt.legend()
    plt.legend(fontsize=20, frameon=True, framealpha=0.8)
    plt.tight_layout()

    if args.out:
        root, ext = os.path.splitext(args.out)
        if ext.lower() == '.pdf':
            plt.savefig(args.out, bbox_inches='tight', transparent=True)
        else:
            plt.savefig(args.out, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {args.out}")
    else:
        plt.show()

if __name__ == '__main__':
    main()
