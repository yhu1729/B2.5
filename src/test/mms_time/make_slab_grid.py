#!/usr/bin/env python
"""Write a rectangular slab grid for b2ag in the pre-2012 free-format
Carre/DG layout (b2ag.dat: parg(0) = -1, structured mode).

Poloidal direction: Z (file index i, nx cells); radial direction: R (file
index j, ny cells). Guard cells are one extra cell ring: file indices run
0..nx+1 and 0..ny+1 (B2 index = file index - 1). Toroidal symmetry
(isymm = 1) at a large major radius keeps the standard code path; the
toroidal field is Bt = R0*Bt0/R, the poloidal field Bp is uniform.

Reader: modules/B2.5/src/preprocessing/b2agfs_st.F, old-format branch
(VERSION < 01.001.028): line 1 'VERSION01.001.000', line 2 'nnx nny', then
one line per cell, iy outer, ix inner:
  i j  xc yc  x1 y1 x2 y2 x3 y3 x4 y4  bp bt
with corners 1..4 = B2 corners 0..3 in the "2-3 / 0-1" convention, corner
0->1 along the poloidal direction.
"""
import argparse


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--nx", type=int, default=48, help="poloidal interior cells")
    ap.add_argument("--ny", type=int, default=12, help="radial interior cells")
    ap.add_argument("--lx", type=float, default=2.0, help="poloidal extent [m] (Z)")
    ap.add_argument("--ly", type=float, default=0.06, help="radial extent [m] (R)")
    ap.add_argument("--r0", type=float, default=10.0, help="major radius of the box centre [m]")
    ap.add_argument("--bt0", type=float, default=2.0, help="toroidal field at r0 [T]")
    ap.add_argument("--bp", type=float, default=0.2, help="uniform poloidal field [T]")
    ap.add_argument("-o", "--output", default="mms_slab.geo")
    a = ap.parse_args()

    dx = a.lx / a.nx
    dy = a.ly / a.ny
    zmax = 0.5 * a.lx
    rmin = a.r0 - 0.5 * a.ly

    # The poloidal index runs towards decreasing Z so that (poloidal, radial)
    # is a right-handed pair with radial = increasing R, as B2 requires
    # (b2agmt_st: qz = sin(angle from x to y direction) must be positive).
    def zc(i):  # poloidal edge coordinate for file index i (cell i spans z(i)..z(i+1), decreasing)
        return zmax - (i - 1) * dx

    def rc(j):
        return rmin + (j - 1) * dy

    fmt = "%.16e"
    lines = ["VERSION01.001.000", "%d %d" % (a.nx, a.ny)]
    for j in range(0, a.ny + 2):
        r_lo, r_hi = rc(j), rc(j + 1)
        r_mid = 0.5 * (r_lo + r_hi)
        for i in range(0, a.nx + 2):
            z_a, z_b = zc(i), zc(i + 1)   # z_a > z_b: corner 0 -> 1 goes towards decreasing Z
            z_mid = 0.5 * (z_a + z_b)
            # B2 corners (2-3 / 0-1): 0 = (r_lo, z_a), 1 = (r_lo, z_b), 2 = (r_hi, z_a), 3 = (r_hi, z_b)
            pts = [(r_mid, z_mid), (r_lo, z_a), (r_lo, z_b), (r_hi, z_a), (r_hi, z_b)]
            bt = a.r0 * a.bt0 / r_mid
            vals = [i, j] + [fmt % v for p in pts for v in p] + [fmt % a.bp, fmt % bt]
            lines.append(" ".join(str(v) for v in vals))
    with open(a.output, "w") as f:
        f.write("\n".join(lines) + "\n")
    print("wrote %s: %d x %d interior cells, dx = %.4g m, dy = %.4g m, R in [%.4g, %.4g], Z in [%.4g, %.4g]"
          % (a.output, a.nx, a.ny, dx, dy, rc(1), rc(a.ny + 1), zc(a.nx + 1), zc(1)))


if __name__ == "__main__":
    main()
