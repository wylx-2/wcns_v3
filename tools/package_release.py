#!/usr/bin/env python3
"""Create or verify a deterministic private WCNS source release archive."""

from __future__ import annotations

import argparse
import gzip
import hashlib
import re
import shutil
import subprocess
import tarfile
import tempfile
from pathlib import Path, PurePosixPath


REPOSITORY = Path(__file__).resolve().parents[1]
VERSION = "1.1.0"

DIRECTORY_PREFIXES = (
    Path("include"),
    Path("src"),
    Path("tests"),
    Path("tools"),
    Path("docs"),
    Path("examples"),
    Path("cases/config"),
    Path("third_party/cgns"),
)

EXPLICIT_FILES = (
    Path("CMakeLists.txt"),
    Path("README.md"),
    Path("LICENSE.md"),
    Path("THIRD_PARTY_NOTICES.md"),
    Path("\u7b97\u6cd5\u8865\u5145.md"),
    Path(".github/workflows/ci.yml"),
)

GENERATED_FILES = (
    PurePosixPath("WCNS_SOURCE_REVISION"),
    PurePosixPath("PACKAGE_CONTENTS.sha256"),
)


def run_git(*arguments: str) -> str:
    completed = subprocess.run(
        ["git", "-c", "core.quotepath=false", *arguments],
        cwd=REPOSITORY,
        text=True,
        encoding="utf-8",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"git {' '.join(arguments)} failed: {completed.stderr.strip()}")
    return completed.stdout.strip()


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--verify",
        type=Path,
        metavar="ARCHIVE",
        help="verify an existing archive and optional adjacent .sha256 sidecar",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=REPOSITORY / "dist",
        help="directory receiving the tar.gz and SHA-256 sidecar",
    )
    return parser.parse_args()


def tracked_payload() -> list[Path]:
    tracked = {
        Path(line)
        for line in run_git("ls-tree", "-r", "--name-only", "HEAD").splitlines()
        if line
    }
    selected = set(EXPLICIT_FILES)
    for path in tracked:
        if any(path == prefix or prefix in path.parents
               for prefix in DIRECTORY_PREFIXES):
            selected.add(path)

    missing = sorted(path for path in selected if path not in tracked)
    if missing:
        raise RuntimeError(
            "release payload contains files not committed in HEAD: "
            + ", ".join(path.as_posix() for path in missing))
    absent = sorted(path for path in selected if not (REPOSITORY / path).is_file())
    if absent:
        raise RuntimeError(
            "release payload is missing from the worktree: "
            + ", ".join(path.as_posix() for path in absent))

    forbidden = [
        path for path in selected
        if path == Path(".git")
        or Path("cases/manual") in path.parents
        or any(part.startswith("build-") for part in path.parts)
    ]
    if forbidden:
        raise RuntimeError(
            "forbidden release payload entries: "
            + ", ".join(path.as_posix() for path in forbidden))

    return sorted(selected, key=lambda item: item.as_posix())


def require_clean_payload(paths: list[Path]) -> None:
    command = ["git", "diff", "--quiet", "HEAD", "--"]
    command.extend(path.as_posix() for path in paths)
    completed = subprocess.run(command, cwd=REPOSITORY, check=False)
    if completed.returncode == 1:
        raise RuntimeError(
            "selected release files differ from HEAD; commit them before packaging")
    if completed.returncode != 0:
        raise RuntimeError("unable to verify release payload against HEAD")


def require_version(cmake_text: str) -> None:
    match = re.search(
        r"project\s*\(\s*wcns\s+VERSION\s+([0-9]+\.[0-9]+\.[0-9]+)",
        cmake_text,
        flags=re.IGNORECASE,
    )
    if match is None or match.group(1) != VERSION:
        observed = "missing" if match is None else match.group(1)
        raise RuntimeError(
            f"CMake project version is {observed}; expected {VERSION}")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_bytes(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def normalized_tar_info(
    info: tarfile.TarInfo, timestamp: int,
) -> tarfile.TarInfo:
    info.uid = 0
    info.gid = 0
    info.uname = ""
    info.gname = ""
    info.mtime = timestamp
    info.pax_headers = {}
    if info.isdir():
        info.mode = 0o755
    elif info.name.endswith(".py"):
        info.mode = 0o755
    else:
        info.mode = 0o644
    return info


def create_archive(
    package_root: Path, archive: Path, timestamp: int,
) -> None:
    entries = [package_root]
    entries.extend(sorted(package_root.rglob("*"), key=lambda path: path.as_posix()))
    with archive.open("wb") as raw:
        with gzip.GzipFile(
            filename="", mode="wb", fileobj=raw, compresslevel=9, mtime=timestamp,
        ) as compressed:
            with tarfile.open(
                fileobj=compressed, mode="w", format=tarfile.PAX_FORMAT,
            ) as tar:
                for entry in entries:
                    relative = entry.relative_to(package_root)
                    arcname = package_root.name
                    if relative.parts:
                        arcname += "/" + relative.as_posix()
                    tar.add(
                        entry,
                        arcname=arcname,
                        recursive=False,
                        filter=lambda info: normalized_tar_info(info, timestamp),
                    )


def parse_contents_manifest(value: bytes) -> dict[PurePosixPath, str]:
    result: dict[PurePosixPath, str] = {}
    try:
        text = value.decode("utf-8")
    except UnicodeDecodeError as error:
        raise RuntimeError("PACKAGE_CONTENTS.sha256 is not UTF-8") from error
    for line_number, line in enumerate(text.splitlines(), start=1):
        match = re.fullmatch(r"([0-9a-f]{64})  (.+)", line)
        if match is None:
            raise RuntimeError(
                f"invalid PACKAGE_CONTENTS.sha256 line {line_number}")
        path = PurePosixPath(match.group(2))
        if path.is_absolute() or ".." in path.parts or path in result:
            raise RuntimeError(
                f"unsafe or duplicate manifest path on line {line_number}")
        result[path] = match.group(1)
    return result


def verify_archive(archive: Path, check_sidecar: bool = True) -> dict[str, object]:
    archive = archive.expanduser().resolve()
    if not archive.is_file():
        raise RuntimeError(f"archive does not exist: {archive}")

    expected_prefix = f"wcns-{VERSION}-source-"
    regular: dict[PurePosixPath, bytes] = {}
    roots: set[str] = set()
    with tarfile.open(archive, mode="r:gz") as tar:
        seen: set[str] = set()
        for member in tar.getmembers():
            member_path = PurePosixPath(member.name)
            if (member_path.is_absolute() or ".." in member_path.parts
                    or member.name in seen):
                raise RuntimeError(f"unsafe or duplicate archive member: {member.name}")
            seen.add(member.name)
            if not member_path.parts:
                raise RuntimeError("archive contains an empty member path")
            roots.add(member_path.parts[0])
            if member.issym() or member.islnk() or member.isdev():
                raise RuntimeError(f"unsupported archive member type: {member.name}")
            if member.isfile():
                stream = tar.extractfile(member)
                if stream is None:
                    raise RuntimeError(f"cannot read archive member: {member.name}")
                relative = PurePosixPath(*member_path.parts[1:])
                regular[relative] = stream.read()

    if len(roots) != 1:
        raise RuntimeError("archive must contain exactly one top-level directory")
    root = next(iter(roots))
    if not root.startswith(expected_prefix) or not re.fullmatch(
            re.escape(expected_prefix) + r"[0-9a-f]{12}", root):
        raise RuntimeError(f"unexpected archive root: {root}")

    for forbidden in regular:
        if (".git" in forbidden.parts
                or PurePosixPath("cases/manual") in forbidden.parents
                or any(part.startswith("build-") for part in forbidden.parts)):
            raise RuntimeError(f"forbidden archive member: {forbidden}")

    if PurePosixPath("PACKAGE_CONTENTS.sha256") not in regular:
        raise RuntimeError("archive lacks PACKAGE_CONTENTS.sha256")
    manifest = parse_contents_manifest(
        regular[PurePosixPath("PACKAGE_CONTENTS.sha256")])
    expected_files = set(regular) - {PurePosixPath("PACKAGE_CONTENTS.sha256")}
    if set(manifest) != expected_files:
        missing = sorted(str(path) for path in expected_files - set(manifest))
        extra = sorted(str(path) for path in set(manifest) - expected_files)
        raise RuntimeError(
            f"content manifest mismatch; missing={missing}, extra={extra}")
    for path, expected in manifest.items():
        observed = sha256_bytes(regular[path])
        if observed != expected:
            raise RuntimeError(f"content hash mismatch: {path}")

    revision = regular.get(PurePosixPath("WCNS_SOURCE_REVISION"), b"").decode(
        "ascii", errors="strict").strip()
    if not re.fullmatch(r"[0-9a-f]{12}", revision):
        raise RuntimeError("invalid WCNS_SOURCE_REVISION")
    if root != f"{expected_prefix}{revision}":
        raise RuntimeError("archive root and WCNS_SOURCE_REVISION disagree")
    cmake = regular.get(PurePosixPath("CMakeLists.txt"))
    if cmake is None:
        raise RuntimeError("archive lacks CMakeLists.txt")
    require_version(cmake.decode("utf-8"))

    archive_hash = sha256(archive)
    sidecar = archive.with_name(archive.name + ".sha256")
    sidecar_checked = False
    if check_sidecar and sidecar.exists():
        expected_sidecar = f"{archive_hash}  {archive.name}\n"
        if sidecar.read_text(encoding="utf-8") != expected_sidecar:
            raise RuntimeError("archive SHA-256 sidecar does not match")
        sidecar_checked = True

    return {
        "archive": str(archive),
        "sha256": archive_hash,
        "revision": revision,
        "payload_files": len(regular),
        "sidecar_checked": sidecar_checked,
    }


def create_release(output_directory: Path) -> tuple[Path, Path]:
    output_directory = output_directory.expanduser().resolve()
    output_directory.mkdir(parents=True, exist_ok=True)

    revision = run_git("rev-parse", "--short=12", "HEAD")
    timestamp = int(run_git("show", "-s", "--format=%ct", "HEAD"))
    payload = tracked_payload()
    require_clean_payload(payload)
    require_version((REPOSITORY / "CMakeLists.txt").read_text(encoding="utf-8"))

    package_name = f"wcns-{VERSION}-source-{revision}"
    archive = output_directory / f"{package_name}.tar.gz"
    sidecar = output_directory / f"{package_name}.tar.gz.sha256"

    with tempfile.TemporaryDirectory(
        prefix="wcns-release-package-", dir=output_directory,
    ) as temporary:
        package_root = Path(temporary) / package_name
        package_root.mkdir()
        for relative in payload:
            destination = package_root / relative
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(REPOSITORY / relative, destination)

        (package_root / "WCNS_SOURCE_REVISION").write_text(
            revision + "\n", encoding="utf-8", newline="\n")
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

        temporary_archive = Path(temporary) / archive.name
        create_archive(package_root, temporary_archive, timestamp)
        shutil.copyfile(temporary_archive, archive)

    archive_hash = sha256(archive)
    sidecar.write_text(
        f"{archive_hash}  {archive.name}\n", encoding="utf-8", newline="\n")
    result = verify_archive(archive)
    print(f"created: {archive}")
    print(f"created: {sidecar}")
    print(f"revision: {revision}")
    print(f"payload files: {result['payload_files']}")
    print(f"archive bytes: {archive.stat().st_size}")
    print(f"sha256: {archive_hash}")
    return archive, sidecar


def main() -> int:
    args = parse_args()
    if args.verify is not None:
        result = verify_archive(args.verify)
        print("verified: " + str(result["archive"]))
        print("revision: " + str(result["revision"]))
        print("payload files: " + str(result["payload_files"]))
        print("sha256: " + str(result["sha256"]))
        print("sidecar checked: " + str(result["sidecar_checked"]).lower())
        return 0
    create_release(args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
