#!/usr/bin/env python3
import argparse, os, csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype']  = 42

def read_cpu_csv(path):
    with open(path, 'r', newline='') as f:
        rdr = csv.DictReader(f)
        hs, lgs = [], []
        cols = {'vanilla':[], 'vanilla_feedback':[], 'feedback':[], 'adaptive':[], 'strang':[]}
        for row in rdr:
            h = float(row['h'])
            hs.append(h)
            lg = float(row['log10_h']) if ('log10_h' in row and row['log10_h'] != '') else np.log10(h)
            lgs.append(lg)
            for k in cols:
                cols[k].append(float(row[k]))
    hs  = np.array(hs)
    lgs = np.array(lgs)
    data = {k: np.array(v) for k, v in cols.items()}
    idx = np.argsort(lgs)
    lgs = lgs[idx]
    for k in data:
        data[k] = data[k][idx]
    return lgs, data

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--csv',  default='results/cpu_time/cpu_time_vs_h_all.csv')
    ap.add_argument('--out',  default=None)
    ap.add_argument('--title', default='CPU time vs Step Size')
    ap.add_argument('--x', choices=['logh','neglogh','N'], default='neglogh',
                    help="x-axis: logh=log10(h), neglogh=-log10(h), N=log10(tf/h)")
    ap.add_argument('--tf', type=float, default=None,
                    help="final time tf; required if --x N")
    args = ap.parse_args()

    x_log10h, d = read_cpu_csv(args.csv)

    if args.x == 'logh':
        x = x_log10h
        xlab = r'$\log_{10} h$'
    elif args.x == 'neglogh':
        x = -x_log10h
        xlab = r'$\log_{10}(1/h)$'
    else:  # 'N'
        if args.tf is None:
            raise SystemExit("ERROR: --tf is required when --x N")
        # h = 10**(x_log10h)  ⇒  N = tf / h = tf * 10**(-x_log10h)
        N = args.tf * (10.0 ** (-x_log10h))
        x = np.log10(N)
        xlab = r'$\log_{10} N\;(N=t_f/h)$'

    fig, ax = plt.subplots(figsize=(12,8), constrained_layout=True)

    order = [
        ('vanilla',          "Euler's method",                'o'),
        ('vanilla_feedback', 'Feedback (unity gain)',         's'),
        ('feedback',         r'Feedback $(1/hL)$',            '^'),
        ('adaptive',         'Adaptive Feedback',             'D'),
        ('SV',               'Splitting method',         'v'),
    ]
    for key, label, marker in order:
        if key in d:
            y = d[key]
            m = np.isfinite(y)
            ax.plot(x[m], y[m], marker=marker, markersize=12,
                    linewidth=2.4, label=label)

    ax.set_xlabel(xlab, fontsize=26)
    ax.set_ylabel('CPU time (s)', fontsize=26)
    ax.grid(True, linestyle='--', alpha=0.8)
    ax.tick_params(labelsize=20)
    ax.legend(fontsize=20, frameon=True, framealpha=0.8)
    ax.set_yscale('log')

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
