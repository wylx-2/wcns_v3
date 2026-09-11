#!/usr/bin/env python3
"""Create a deterministic, minimal Linux source package for WCNS case05."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path


REPOSITORY = Path(__file__).resolve().parents[1]
CASE = Path("cases/manual/case05_3d_turbulent_channel")

DIRECTORY_PREFIXES = (
    Path("include"),
    Path("src"),
)

EXPLICIT_FILES = (
    Path("CMakeLists.txt"),
    Path("LICENSE.md"),
    Path("THIRD_PARTY_NOTICES.md"),
    Path("third_party/cgns/CGNS-4.4.0.zip"),
    Path("third_party/cgns/README.md"),
    Path("tools/generate_release_cgns.cpp"),
    Path("tools/validate_release_case.cpp"),
    Path("tools/compare_metric_profiles.cpp"),
    CASE / "LINUX_SERVER_GUIDE.md",
    CASE / "channel_retau180_feasibility.wcns",
    CASE / "channel_retau180_longrun.wcns",
    CASE / "run_case05.py",
    CASE / "validate_case05.py",
)


def run_git(*arguments: str) -> str:
    completed = subprocess.run(
        ["git", *arguments],
        cwd=REPOSITORY,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"git {' '.join(arguments)} failed: {completed.stderr.strip()}")
    return completed.stdout.strip()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY / "dist",
        help="directory receiving the tar.gz and its SHA-256 sidecar",
    )
    return parser.parse_args()


def tracked_payload() -> list[Path]:
    tracked = {
        Path(line) for line in run_git("ls-tree", "-r", "--name-only", "HEAD").splitlines() if line
    }
    selected = set(EXPLICIT_FILES)
    for path in tracked:
        if any(path == prefix or prefix in path.parents for prefix in DIRECTORY_PREFIXES):
            selected.add(path)

    missing = sorted(path for path in selected if path not in tracked)
    if missing:
        raise RuntimeError(
            "package payload contains files not committed in HEAD: "
            + ", ".join(path.as_posix() for path in missing)
        )
    absent = sorted(path for path in selected if not (REPOSITORY / path).is_file())
    if absent:
        raise RuntimeError(
            "package payload is missing from the worktree: "
            + ", ".join(path.as_posix() for path in absent)
        )
    return sorted(selected, key=lambda item: item.as_posix())


def require_clean_payload(paths: list[Path]) -> None:
    command = ["git", "diff", "--quiet", "HEAD", "--"]
    command.extend(path.as_posix() for path in paths)
    completed = subprocess.run(command, cwd=REPOSITORY, check=False)
    if completed.returncode == 1:
        raise RuntimeError("selected package files differ from HEAD; commit them before packaging")
    if completed.returncode != 0:
        raise RuntimeError("unable to verify package payload against HEAD")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_tar_info(
    info: tarfile.TarInfo,
    timestamp: int,
) -> tarfile.TarInfo:
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    info.mtime = timestamp
    if info.isdir():
        info.mode = 0o755
    elif info.name.endswith(".py"):
        info.mode = 0o755
    else:
        info.mode = 0o644
    return info


def create_archive(source: Path, archive: Path, timestamp: int) -> None:
    with archive.open("wb") as raw:
        with gzip.GzipFile(
            filename="",
            mode="wb",
            fileobj=raw,
            compresslevel=9,
            mtime=timestamp,
        ) as compressed:
            with tarfile.open(
                fileobj=compressed,
                mode="w",
                format=tarfile.PAX_FORMAT,
            ) as tar:
                tar.add(
                    source,
                    arcname=source.name,
                    recursive=True,
                    filter=lambda info: normalized_tar_info(info, timestamp),
                )


def main() -> int:
    args = parse_args()
    output_directory = args.output_dir.expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    revision = run_git("rev-parse", "--short=12", "HEAD")
    timestamp = int(run_git("show", "-s", "--format=%ct", "HEAD"))
    payload = tracked_payload()
    require_clean_payload(payload)

    package_name = f"wcns-case05-linux-{revision}"
    archive = output_directory / f"{package_name}.tar.gz"
    sidecar = output_directory / f"{package_name}.tar.gz.sha256"

    with tempfile.TemporaryDirectory(
        prefix="wcns-case05-package-",
        dir=output_directory,
    ) as temporary:
        package_root = Path(temporary) / package_name
        package_root.mkdir()
        for relative in payload:
            destination = package_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(REPOSITORY / relative, destination)

        (package_root / "WCNS_SOURCE_REVISION").write_text(
            revision + "\n", encoding="utf-8", newline="\n"
        )

        manifest_paths = [Path("WCNS_SOURCE_REVISION"), *payload]
        manifest_lines = [
            f"{sha256(package_root / path)}  {path.as_posix()}"
            for path in sorted(manifest_paths, key=lambda item: item.as_posix())
        ]
        (package_root / "PACKAGE_CONTENTS.sha256").write_text(
            "\n".join(manifest_lines) + "\n",
            encoding="utf-8",
            newline="\n",
        )
        create_archive(package_root, archive, timestamp)

    archive_hash = sha256(archive)
    sidecar.write_text(f"{archive_hash}  {archive.name}\n", encoding="utf-8", newline="\n")
    print(f"created: {archive}")
    print(f"created: {sidecar}")
    print(f"revision: {revision}")
    print(f"payload files: {len(payload) + 2}")
    print(f"archive bytes: {archive.stat().st_size}")
    print(f"sha256: {archive_hash}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
