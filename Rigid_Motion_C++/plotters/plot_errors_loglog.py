#!/usr/bin/env python3
import argparse, os, csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42

METHOD_LABELS = [
    ("vanilla",            "Euler's method",             'o'),
    ("vanilla_feedback",   "Feedback (unity gain)",      's'),
    ("feedback",           r"Feedback $(1/hL)$",         '^'),
    ("adaptive",           "Adaptive Feedback",          'D'),
    ("strang",             "Splitting method",           'v'),
]

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

def positive_mask(*arrays):
    """Return mask where all arrays are strictly positive at the same indices."""
    mask = np.ones_like(arrays[0], dtype=bool)
    for a in arrays:
        mask &= (a > 0)
    return mask

def plot_panel(ax, h, series_list, colors, ylabel, title, drop_over):
    """
    series_list: [(y_array, key, label, marker), ...]
    colors: dict to keep method->color consistent across panels
    """
    handles = []
    for (y, key, label, marker) in series_list:
        y = np.asarray(y, dtype=float)
        y = np.where(y > drop_over, np.inf, y)
        m = np.isfinite(y) & (y > 0) & np.isfinite(h) & (h > 0)
        if not np.any(m):
            continue

        if key not in colors:
            ln, = ax.loglog(h[m], y[m], marker=marker, markersize=9, linewidth=2.2, label=label)
            colors[key] = ln.get_color()
            handles.append(ln)
        else:
            ln, = ax.loglog(h[m], y[m], marker=marker, markersize=9, linewidth=2.2,
                            label=label, color=colors[key])
            handles.append(ln)

    ax.set_xlabel(r"$h$", fontsize=15)
    ax.set_ylabel(ylabel, fontsize=15)
    if title:
        ax.set_title(title, fontsize=15)
    ax.grid(True, which='both', linestyle='--', alpha=0.8)
    ax.tick_params(labelsize=11)
    return handles


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--errors_csv',  default='results/error_data/errors_vs_h.csv',
                    help='CSV for e (existing error)')
    ap.add_argument('--maxdE_csv',   default='results/error_data/maxdE_vs_h.csv',
                    help='CSV for maxdE')
    ap.add_argument('--maxdPi_csv',  default='results/error_data/maxdPi_vs_h.csv',
                    help='CSV for maxdPi')
    ap.add_argument('--maxdDet_csv', default='results/error_data/maxdDet_vs_h.csv',
                    help='CSV for maxdDet')
    ap.add_argument('--out',   default=None,
                    help='Save figure to this path (.png/.pdf). If omitted, show on screen.')
    ap.add_argument('--title', default='Error & Invariant Deviations vs Step Size',
                    help='Suptitle for the 2x2 figure')
    ap.add_argument('--drop_over', type=float, default=1e8,
                help='Treat values above this as +inf and drop them on log axes')
    args = ap.parse_args()

    # ---- Read all CSVs
    h_err, e_van, e_vanfb, e_fb, e_ad, e_str = read_wide_csv(args.errors_csv)
    h_dE,  dE_van, dE_vanfb, dE_fb, dE_ad, dE_str = read_wide_csv(args.maxdE_csv)
    h_dP,  dP_van, dP_vanfb, dP_fb, dP_ad, dP_str = read_wide_csv(args.maxdPi_csv)
    h_dD,  dD_van, dD_vanfb, dD_fb, dD_ad, dD_str = read_wide_csv(args.maxdDet_csv)

    # ---- Prepare 2x2 figure
    fig, axs = plt.subplots(2, 2, figsize=(15, 10), constrained_layout=True)
    colors = {}  # key -> color, fixed by the first panel and reused

    # Panel (a)
    series_err = [
        (e_van,   "vanilla",           METHOD_LABELS[0][1], METHOD_LABELS[0][2]),
        (e_vanfb, "vanilla_feedback",  METHOD_LABELS[1][1], METHOD_LABELS[1][2]),
        (e_fb,    "feedback",          METHOD_LABELS[2][1], METHOD_LABELS[2][2]),
        (e_ad,    "adaptive",          METHOD_LABELS[3][1], METHOD_LABELS[3][2]),
        (e_str,   "strang",            METHOD_LABELS[4][1], METHOD_LABELS[4][2]),
    ]
    h_a = h_err
    handles = plot_panel(axs[0,0], h_a, series_err, colors,
                         title=None, ylabel=r"$\max \ V(x_k)$", drop_over=args.drop_over)

    # Panel (b)
    series_dD = [
        (dD_van,   "vanilla",          METHOD_LABELS[0][1], METHOD_LABELS[0][2]),
        (dD_vanfb, "vanilla_feedback", METHOD_LABELS[1][1], METHOD_LABELS[1][2]),
        (dD_fb,    "feedback",         METHOD_LABELS[2][1], METHOD_LABELS[2][2]),
        (dD_ad,    "adaptive",         METHOD_LABELS[3][1], METHOD_LABELS[3][2]),
        (dD_str,   "strang",           METHOD_LABELS[4][1], METHOD_LABELS[4][2]),
    ]
    plot_panel(axs[0,1], h_dD, series_dD, colors,
               title=None, ylabel=r"$\max \ \|R^\top R - 1\|_F$", drop_over=args.drop_over)

    # Panel (c)
    series_dE = [
        (dE_van,   "vanilla",          METHOD_LABELS[0][1], METHOD_LABELS[0][2]),
        (dE_vanfb, "vanilla_feedback", METHOD_LABELS[1][1], METHOD_LABELS[1][2]),
        (dE_fb,    "feedback",         METHOD_LABELS[2][1], METHOD_LABELS[2][2]),
        (dE_ad,    "adaptive",         METHOD_LABELS[3][1], METHOD_LABELS[3][2]),
        (dE_str,   "strang",           METHOD_LABELS[4][1], METHOD_LABELS[4][2]),
    ]
    plot_panel(axs[1,0], h_dE, series_dE, colors,
               title=None, ylabel=r"$\max \ |\Delta E|$", drop_over=args.drop_over)

    # Panel (d)
    series_dP = [
        (dP_van,   "vanilla",          METHOD_LABELS[0][1], METHOD_LABELS[0][2]),
        (dP_vanfb, "vanilla_feedback", METHOD_LABELS[1][1], METHOD_LABELS[1][2]),
        (dP_fb,    "feedback",         METHOD_LABELS[2][1], METHOD_LABELS[2][2]),
        (dP_ad,    "adaptive",         METHOD_LABELS[3][1], METHOD_LABELS[3][2]),
        (dP_str,   "strang",           METHOD_LABELS[4][1], METHOD_LABELS[4][2]),
    ]
    plot_panel(axs[1,1], h_dP, series_dP, colors,
               title=None, ylabel=r"$\max \ |\Delta \pi|$", drop_over=args.drop_over)


    # One shared legend (methods), using handles from panel (a)
    if handles:
        fig.legend(handles, [hl.get_label() for hl in handles],
                   loc='outside lower center', ncol=3, frameon=True, framealpha=0.8, fontsize=16)

    # Suptitle
    fig.suptitle(args.title, fontsize=20)

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
