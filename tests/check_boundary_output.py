import glob
import math
import os
import sys


def rows(path):
    with open(path, encoding="utf-8") as stream:
        return [
            [float(value) for value in line.split()]
            for line in stream
            if line.strip() and not line.startswith("#")
        ]


output = sys.argv[1]
rank = int(sys.argv[2]) if len(sys.argv) > 2 else 1
faces = sorted(glob.glob(os.path.join(output, f"*.boundary.r{rank}.step*.txt")))
loads = glob.glob(os.path.join(output, f"*.loads.r{rank}.txt"))
if len(faces) < 2 or len(loads) != 1:
    raise SystemExit(
        f"expected boundary snapshots and one load history, got {faces}, {loads}"
    )

initial_paths = [path for path in faces if ".step00000000." in path]
if len(initial_paths) != 1:
    raise SystemExit(f"expected one initial boundary snapshot, got {initial_paths}")
initial = rows(initial_paths[0])
if not initial:
    raise SystemExit("initial boundary snapshot is empty")
keys = [tuple(row[:7]) for row in initial]
ordered_keys = sorted(
    keys, key=lambda key: (key[0], key[1], key[4], key[3], key[2], key[5], key[6])
)
if keys != ordered_keys or len(keys) != len(set(keys)):
    raise SystemExit("boundary face keys are not sorted and unique")

# Fixed columns end at nz=13; requested Cp is the fourth selected quantity.
cp = [row[17] for row in initial]
if max(abs(value) for value in cp) > 1.0e-13:
    raise SystemExit(f"uniform initial Cp is not zero: {max(abs(v) for v in cp)}")
for row in initial:
    # Requested pressure traction x/y and total traction x/y are columns 18..21.
    if abs(row[18] - row[20]) > 1.0e-13 or abs(row[19] - row[21]) > 1.0e-13:
        raise SystemExit("inviscid face pressure and total traction differ")

history = rows(loads[0])
if len(history) != 2:
    raise SystemExit(f"expected two load rows, got {len(history)}")
first = history[0]
if len(first) != 29 or not all(math.isfinite(value) for value in first):
    raise SystemExit("load history schema or finiteness check failed")
pressure_force = first[2:5]
if math.sqrt(sum(value * value for value in pressure_force)) > 1.0e-12:
    raise SystemExit(f"closed constant-pressure force is nonzero: {pressure_force}")
if any(first[index] != 0.0 for index in range(5, 8)):
    raise SystemExit("inviscid viscous-force columns are not exact zero")
for component in range(3):
    if abs(first[2 + component] + first[5 + component] - first[8 + component]) > 1.0e-13:
        raise SystemExit("pressure plus viscous load does not equal total load")

print("boundary output validation passed")
