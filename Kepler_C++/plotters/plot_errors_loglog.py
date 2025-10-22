#!/usr/bin/env python3
import argparse, os, csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42

METHODS = [
    ("vanilla",           "Euler's method",            'o'),
    ("vanilla_feedback",  "Feedback (unity gain)",     's'),
    ("feedback",          r"Feedback $(1/hL)$",        '^'),
    ("adaptive",          "Adaptive Feedback",         'D'),
    ("SV",                "Störmer–Verlet (B)",        'v'),
]

def read_wide_csv(path):
    with open(path, 'r', newline='') as f:
        rdr = csv.DictReader(f)
        hs, v = [], {k: [] for k,_,_ in METHODS}
        for row in rdr:
            hs.append(float(row['h']))
            for k,_,_ in METHODS:
                v[k].append(float(row[k]))
    return np.array(hs), {k: np.array(v[k]) for k,_,_ in METHODS}

def pos_mask(*arrs):
    m = np.ones_like(arrs[0], dtype=bool)
    for a in arrs:
        m &= (a > 0)
    return m

def plot_panel(ax, h, values_by_key, colors, ylabel, title):
    handles = []
    for key, label, marker in METHODS:
        y = values_by_key[key]
        m = pos_mask(h, y)
        if not np.any(m):  # nothing to draw for this series
            continue
        if key not in colors:
            ln, = ax.loglog(h[m], y[m], marker=marker, markersize=10, linewidth=2.4, label=label)
            colors[key] = ln.get_color()
            handles.append(ln)
        else:
            ln, = ax.loglog(h[m], y[m], marker=marker, markersize=10, linewidth=2.4,
                            label=label, color=colors[key])
            handles.append(ln)
    ax.set_xlabel(r"$h$", fontsize=20)
    ax.set_ylabel(ylabel, fontsize=20)
    ax.set_title(title, fontsize=20)
    ax.grid(True, which='both', linestyle='--', alpha=0.8)
    ax.tick_params(labelsize=16)
    return handles

def mpl_ge_38():
    try:
        major, minor = (int(x) for x in mpl.__version__.split('.')[:2])
        return (major, minor) >= (3, 8)
    except Exception:
        return False

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--errors_csv', default='results/error_data/errors_vs_h.csv')
    ap.add_argument('--maxdL_csv',  default='results/error_data/maxdL_vs_h.csv')
    ap.add_argument('--maxdA_csv',  default='results/error_data/maxdA_vs_h.csv')
    ap.add_argument('--out',        default=None, help='Output path (.png/.pdf). If omitted, shows interactively.')
    ap.add_argument('--title',      default='Error, max dL, max dA vs step size')
    args = ap.parse_args()

    # Read CSVs
    h_err,  err_by_key  = read_wide_csv(args.errors_csv)
    h_dL,   dL_by_key   = read_wide_csv(args.maxdL_csv)
    h_dA,   dA_by_key   = read_wide_csv(args.maxdA_csv)

    # Figure with 1x3 panels
    fig, axs = plt.subplots(1, 3, figsize=(18, 6), constrained_layout=True)
    colors = {}

    # (1) error
    handles = plot_panel(axs[0], h_err, err_by_key, colors,
                         ylabel=r"$\max \ V(x_k)$", title=None)

    # (2) max dL
    plot_panel(axs[1], h_dL, dL_by_key, colors,
               ylabel=r"$\max \ |\Delta L|$", title=None)

    # (3) max dA
    plot_panel(axs[2], h_dA, dA_by_key, colors,
               ylabel=r"$\max \ |\Delta A|$", title=None)

    # Shared legend outside (avoid overlap with axes; no tight_layout/subplots_adjust)
    labels = [h.get_label() for h in handles] if handles else []
    if mpl_ge_38():
        fig.legend(handles, labels, loc='outside lower center', ncol=3,
                   frameon=True, framealpha=0.8, fontsize=16)
    else:
        fig.legend(handles, labels, loc='lower center', ncol=3,
                   bbox_to_anchor=(0.5, -0.03), bbox_transform=fig.transFigure,
                   frameon=True, framealpha=0.8, fontsize=16)

    fig.suptitle(args.title, fontsize=20)

    if args.out:
        root, ext = os.path.splitext(args.out)
        if ext.lower() == '.pdf':
            fig.savefig(args.out, bbox_inches='tight', transparent=True)
        else:
            fig.savefig(args.out, dpi=300, bbox_inches='tight')
        print(f"Saved figure to {args.out}")
    else:
        plt.show()

if __name__ == '__main__':
    main()
