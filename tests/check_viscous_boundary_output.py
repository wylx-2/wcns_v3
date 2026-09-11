import glob
import os
import sys


def relative_error(actual, expected):
    return abs(actual - expected) / max(abs(expected), 1.0e-300)


case = sys.argv[1]
directory = sys.argv[2]
rank = int(sys.argv[3])
paths = glob.glob(os.path.join(directory, f"*.boundary.r{rank}.step00000000.*.txt"))
if len(paths) != 1:
    raise SystemExit(f"expected one initial boundary snapshot, got {paths}")
with open(paths[0], encoding="utf-8") as stream:
    rows = [
        [float(value) for value in line.split()]
        for line in stream
        if line.strip() and not line.startswith("#")
    ]
if not rows:
    raise SystemExit("viscous boundary snapshot is empty")

lower = [row for row in rows if int(row[0]) in (0, 1)]
upper = [row for row in rows if int(row[0]) in (2, 3)]
if len(lower) + len(upper) != len(rows) or not lower or not upper:
    raise SystemExit("viscous boundary patch coverage is incomplete")

# Base columns: 0..13. Selected: p,T,mu,Cp,Cf,q,viscous_tx,...
reynolds = 1.225 * 340.0 / 1.7894e-5
mach2 = 340.0**2 / (1.4 * (8.314 / 0.029) * 288.15)
chi = 1.0 / ((1.4 - 1.0) * mach2 * 0.72)
if case == "couette":
    expected_cf = (0.2 / reynolds) / 0.5
    for row in lower:
        if relative_error(row[18], expected_cf) > 5.0e-11:
            raise SystemExit(f"lower Couette Cf mismatch: {row[18]} vs {expected_cf}")
    for row in upper:
        if relative_error(row[18], -expected_cf) > 5.0e-11:
            raise SystemExit(f"upper Couette Cf mismatch: {row[18]} vs {-expected_cf}")
    if max(abs(row[19]) for row in rows) > 1.0e-13:
        raise SystemExit("isothermal constant-temperature Couette heat flux is nonzero")
elif case == "conduction":
    expected_heat = chi / reynolds
    for row in lower:
        if relative_error(row[19], expected_heat) > 5.0e-11:
            raise SystemExit(f"lower conduction heat flux mismatch: {row[19]} vs {expected_heat}")
    for row in upper:
        if relative_error(row[19], -expected_heat) > 5.0e-11:
            raise SystemExit(f"upper conduction heat flux mismatch: {row[19]} vs {-expected_heat}")
    if max(abs(row[18]) for row in rows) > 1.0e-13:
        raise SystemExit("linear conduction skin friction is nonzero")
else:
    raise SystemExit(f"unknown viscous boundary case: {case}")

print(f"{case} boundary trace validation passed for rank {rank}")
