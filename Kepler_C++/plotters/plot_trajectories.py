#!/usr/bin/env python3
import argparse, os, csv
import numpy as np
import matplotlib.pyplot as plt

def parse_range(s):
    a, b = s.split(',')
    return float(a), float(b)

def find_csv(base, method, h_str):
    folder = os.path.join(base, method, f"h={h_str}")
    csv_path = os.path.join(folder, "traj.csv")
    if os.path.exists(csv_path):
        return csv_path

    # fallback: 가장 가까운 h 디렉터리 선택 (log-스케일 거리)
    mdir = os.path.join(base, method)
    if not os.path.isdir(mdir):
        raise SystemExit(f"No method directory: {mdir}")
    candidates = [d for d in os.listdir(mdir) if d.startswith("h=")]
    if not candidates:
        raise SystemExit(f"No trajectories found in {mdir}")
    target = float(h_str)
    def closeness(d):
        val = float(d.split('=')[1])
        return abs(np.log10(val) - np.log10(target))
    folder = os.path.join(mdir, sorted(candidates, key=closeness)[0])
    return os.path.join(folder, "traj.csv")

def load_x1x2(csv_path):
    with open(csv_path, 'r') as f:
        rdr = csv.DictReader(f)
        fieldnames = [c.strip() for c in rdr.fieldnames] if rdr.fieldnames else []
        has_x1x2 = ('x1' in fieldnames and 'x2' in fieldnames)

        x1, x2 = [], []
        if has_x1x2:
            for row in rdr:
                x1.append(float(row['x1']))
                x2.append(float(row['x2']))
            return np.array(x1), np.array(x2)
        else:
            raise SystemExit(
                f"Required columns not found. Got columns: {fieldnames}\n"
            )

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', default='results/trajectory_data')
    ap.add_argument('--h', required=True, help='e.g., 0.001 or 1e-3 (string match in folder name)')
    ap.add_argument('--method', default='vanilla',
                    choices=['vanilla','vanilla_feedback','feedback','adaptive','SV'])
    ap.add_argument('--out', default=None, help='save path (PDF/PNG). If omitted, show interactively.')
    ap.add_argument('--title', default=None)
    ap.add_argument('--xlim', default='auto', help='xmin,xmax or "auto"')
    ap.add_argument('--ylim', default='auto', help='ymin,ymax or "auto"')
    ap.add_argument('--mark_endpoints', type=int, default=False, help='plot start/end markers (1/0)')
    ap.add_argument('--linewidth', type=float, default=3.0)
    args = ap.parse_args()

    csv_path = find_csv(args.base, args.method, args.h)
    x1, x2 = load_x1x2(csv_path)

    fig = plt.figure(figsize=(12.0, 8.0), constrained_layout=True)
    ax = fig.add_subplot(111)
    ax.plot(x1, x2, linewidth=args.linewidth, color='k')
    if args.mark_endpoints:
        ax.plot(x1[0], x2[0], marker='o', markersize=6)   # start
        ax.plot(x1[-1], x2[-1], marker='x', markersize=8) # end

    ax.set_xlabel(r'$x_1$', fontsize=40)
    ax.set_ylabel(r'$x_2$', fontsize=40)
    ax.tick_params(axis='both', labelsize=36)
    ax.grid(True, which='both', linestyle='--', alpha=0.8)
    # ax.set_aspect('equal', adjustable='box')
    if args.title:
        ax.set_title(args.title, fontsize=36)

    # axis limits
    if args.xlim != 'auto':
        ax.set_xlim(*parse_range(args.xlim))
    if args.ylim != 'auto':
        ax.set_ylim(*parse_range(args.ylim))

    if args.out:
        os.makedirs(os.path.dirname(args.out), exist_ok=True)
        plt.savefig(args.out, bbox_inches='tight', transparent=True)
        print(f"Saved x1-x2 trajectory to {args.out}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
