#!/usr/bin/env python3
import argparse, os, csv
import numpy as np
import matplotlib.pyplot as plt

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
    # equal aspect for 3D
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

def parse_range(s):
    a,b = s.split(',')
    return float(a), float(b)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--base', default='trajectory_data')
    ap.add_argument('--h', required=True, help='e.g., 0.001 or 1e-3 (string match in folder name)')
    ap.add_argument('--method', default='vanilla',
                    choices=['vanilla','vanilla_feedback','feedback','adaptive','strang'])
    ap.add_argument('--out', default=None, help='optional save path (PNG). If omitted, just show.')
    ap.add_argument('--title', default=None, help='optional plot title')
    ap.add_argument('--elev', type=float, default=25.0, help='elevation angle (deg)')
    ap.add_argument('--azim', type=float, default=-125.0, help='azimuth angle (deg)')
    ap.add_argument('--xlim', default='0.5,1.5', help='xmin,xmax or "auto"')
    ap.add_argument('--ylim', default='-2.0,2.0', help='ymin,ymax or "auto"')
    ap.add_argument('--zlim', default='-2.0,2.0', help='zmin,zmax or "auto"')
    args = ap.parse_args()

    folder = os.path.join(args.base, args.method, f"h={args.h}")
    csv_path = os.path.join(folder, "traj.csv")
    if not os.path.exists(csv_path):
        # try approximate: pick first directory starting with f"h="
        mdir = os.path.join(args.base, args.method)
        candidates = [d for d in os.listdir(mdir) if d.startswith("h=")]
        if not candidates:
            raise SystemExit(f"No trajectories found in {mdir}")
        # choose the closest in log-space
        target = float(args.h)
        def closeness(d):
            val = float(d.split('=')[1])
            return abs(np.log10(val) - np.log10(target))
        folder = os.path.join(mdir, sorted(candidates, key=closeness)[0])
        csv_path = os.path.join(folder, "traj.csv")

    t, om, R = load_traj(csv_path)

    # Plot Omega trajectory
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot(om[:,0], om[:,1], om[:,2])
    ax.set_xlabel(r'$\Omega_1$')
    ax.set_ylabel(r'$\Omega_2$')
    ax.set_zlabel(r'$\Omega_3$')
    if args.title:
        ax.set_title(args.title)
    # set_equal_3d(ax, om[:,0], om[:,1], om[:,2])
    ax.view_init(elev=args.elev, azim=args.azim)

    # axis limits: use fixed defaults (≈ the figure) unless "auto"
    any_auto = False
    if args.xlim == 'auto': any_auto = True
    else: ax.set_xlim(*parse_range(args.xlim))
    if args.ylim == 'auto': any_auto = True
    else: ax.set_ylim(*parse_range(args.ylim))
    if args.zlim == 'auto': any_auto = True
    else: ax.set_zlim(*parse_range(args.zlim))

    # if any axis is 'auto', keep equal aspect based on data
    if any_auto:
        set_equal_3d(ax, om[:,0], om[:,1], om[:,2])

    if args.out:
        plt.savefig(args.out, dpi=300, bbox_inches='tight')
        print(f"Saved 3D trajectory to {args.out}")
    else:
        plt.show()

    # # Plot Omega
    # plt.figure()
    # plt.plot(t, om[:,0], label='omega_x')
    # plt.plot(t, om[:,1], label='omega_y')
    # plt.plot(t, om[:,2], label='omega_z')
    # plt.xlabel('t')
    # plt.ylabel('Omega components')
    # plt.legend()
    # if args.out:
    #     base,ext = os.path.splitext(args.out)
    #     out1 = base+"_omega.png"
    #     plt.savefig(out1, dpi=150, bbox_inches='tight')
    # else:
    #     plt.show()

    # # Plot R entries (9 curves)
    # plt.figure()
    # labels = ['R00','R01','R02','R10','R11','R12','R20','R21','R22']
    # for i in range(9):
    #     plt.plot(t, R[:,i], label=labels[i])
    # plt.xlabel('t')
    # plt.ylabel('R entries')
    # plt.legend(ncol=3)
    # if args.out:
    #     base,ext = os.path.splitext(args.out)
    #     out2 = base+"_R.png"
    #     plt.savefig(out2, dpi=150, bbox_inches='tight')
    #     print(f"Saved plots to {out1} and {out2}")
    # else:
    #     plt.show()

if __name__ == "__main__":
    main()
