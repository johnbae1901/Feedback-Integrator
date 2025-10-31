#!/usr/bin/env python3
import argparse, os, csv
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox

def load_traj(csv_path):
    t=[]; om=[]; R=[]
    with open(csv_path, 'r') as f:
        rdr = csv.DictReader(f)
        for row in rdr:
            t.append(float(row['t']))
            om.append([float(row['omega_x']), float(row['omega_y']), float(row['omega_z'])])
            R.append([float(row['R00']),float(row['R01']),float(row['R02']),
                      float(row['R10']),float(row['R11']),float(row['R12']),
                      float(row['R20']),float(row['R21']),float(row['R22'])])
    return np.array(t), np.array(om), np.array(R)

def set_equal_3d(ax, X, Y, Z):
    x_min, x_max = np.min(X), np.max(X)
    y_min, y_max = np.min(Y), np.max(Y)
    z_min, z_max = np.min(Z), np.max(Z)
    x_mid = 0.5*(x_min+x_max); y_mid = 0.5*(y_min+y_max); z_mid = 0.5*(z_min+z_max)
    max_range = max(x_max-x_min, y_max-y_min, z_max-z_min)
    r = 0.5*max_range if max_range > 0 else 1.0
    ax.set_xlim(x_mid-r, x_mid+r)
    ax.set_ylim(y_mid-r, y_mid+r)
    ax.set_zlim(z_mid-r, z_mid+r)
    try:
        ax.set_box_aspect((1,1,1))
    except Exception:
        pass

def save_with_asymmetric_margins(fig, ax, outpath,
                                 pad_inch=(0.18, 0.06, 0.18, 0.06)):
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bb = ax.get_tightbbox(renderer).transformed(fig.dpi_scale_trans.inverted())

    L, B, R, T = pad_inch
    crop = Bbox.from_extents(bb.x0 - L, bb.y0 - B, bb.x1 + R, bb.y1 + T)

    fig.savefig(outpath,
                bbox_inches=crop,
                dpi=300,
                transparent=True)

def parse_range(s):
    a,b = s.split(',')
    return float(a), float(b)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', default='trajectory_data')
    ap.add_argument('--h', required=True)
    ap.add_argument('--method', default='vanilla',
                    choices=['vanilla','vanilla_feedback','feedback','adaptive','strang'])
    ap.add_argument('--out', default=None, help='PNG path')
    ap.add_argument('--title', default=None)
    ap.add_argument('--elev', type=float, default=25.0)
    ap.add_argument('--azim', type=float, default=-130.0)
    ap.add_argument('--xlim', default='0.5,1.5')   # or "auto"
    ap.add_argument('--ylim', default='-2.0,2.0')  # or "auto"
    ap.add_argument('--zlim', default='-2.0,2.0')  # or "auto"
    args = ap.parse_args()

    folder = os.path.join(args.base, args.method, f"h={args.h}")
    csv_path = os.path.join(folder, "traj.csv")
    if not os.path.exists(csv_path):
        mdir = os.path.join(args.base, args.method)
        candidates = [d for d in os.listdir(mdir) if d.startswith("h=")]
        if not candidates:
            raise SystemExit(f"No trajectories found in {mdir}")
        target = float(args.h)
        def closeness(d):
            val = float(d.split('=')[1])
            return abs(np.log10(val) - np.log10(target))
        folder = os.path.join(mdir, sorted(candidates, key=closeness)[0])
        csv_path = os.path.join(folder, "traj.csv")

    t, om, R = load_traj(csv_path)

    # === Plot Omega trajectory (3D) ===
    fig = plt.figure(figsize=(8.0,8.0), constrained_layout=True)
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(om[:,0], om[:,1], om[:,2], linewidth=3.0, color='k')

    # labels — add extra padding so they don't collide with tick labels
    ax.xaxis.set_rotate_label(False)
    ax.yaxis.set_rotate_label(False)
    ax.zaxis.set_rotate_label(False)
    ax.set_xlabel(r'$\Omega_1$', fontsize=26, labelpad=12)
    ax.set_ylabel(r'$\Omega_2$', fontsize=26, labelpad=10)
    ax.set_zlabel(r'$\Omega_3$', rotation=0, fontsize=26, labelpad=10)

    # tick label size + distance from axes
    ax.tick_params(axis='x', which='major', labelsize=18, pad=3)
    ax.tick_params(axis='y', which='major', labelsize=18, pad=3)
    ax.tick_params(axis='z', which='major', labelsize=18, pad=3)

    if args.title:
        ax.set_title(args.title)

    ax.view_init(elev=args.elev, azim=args.azim)

    # axis limits
    any_auto = False
    if args.xlim == 'auto': any_auto = True
    else: ax.set_xlim(*parse_range(args.xlim))
    if args.ylim == 'auto': any_auto = True
    else: ax.set_ylim(*parse_range(args.ylim))
    if args.zlim == 'auto': any_auto = True
    else: ax.set_zlim(*parse_range(args.zlim))
    if any_auto:
        set_equal_3d(ax, om[:,0], om[:,1], om[:,2])

    # minimal ticks as requested
    ax.set_xticks([0.5, 1.0, 1.5])
    ax.set_yticks([-2.0, -1.0, 0.0, 1.0, 2.0])
    ax.set_zticks([-2.0, -1.0, 0.0, 1.0, 2.0])

    ax.set_proj_type('ortho')          

    ax.tick_params(axis='x', pad=3)
    ax.tick_params(axis='y', pad=15)
    ax.tick_params(axis='z', pad=3)

    for t in ax.get_xticklabels():
        t.set_verticalalignment('top')
        t.set_horizontalalignment('center')

    for t in ax.get_yticklabels():
        t.set_verticalalignment('bottom') 
        t.set_horizontalalignment('right')

    if args.out:
        fig.canvas.draw()  # update tight bbox
        save_with_asymmetric_margins(fig, ax, args.out,
                             pad_inch=(0.6, 0.0, 0.3, 0.0))
        print(f"Saved 3D trajectory to {args.out}")
    else:
        plt.show()

if __name__ == "__main__":
    main()
