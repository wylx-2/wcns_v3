#!/usr/bin/env python3
"""Small dependency-free process-tree RSS sampler for release measurements."""

from __future__ import annotations

import os
import subprocess
import time
from pathlib import Path


def _windows_process_table() -> dict[int, tuple[int, str]]:
    import ctypes
    from ctypes import wintypes

    class ProcessEntry(ctypes.Structure):
        _fields_ = [
            ("size", wintypes.DWORD),
            ("usage", wintypes.DWORD),
            ("process_id", wintypes.DWORD),
            ("default_heap", ctypes.c_void_p),
            ("module_id", wintypes.DWORD),
            ("threads", wintypes.DWORD),
            ("parent_process_id", wintypes.DWORD),
            ("base_priority", wintypes.LONG),
            ("flags", wintypes.DWORD),
            ("executable", wintypes.WCHAR * 260),
        ]

    snapshot_flag = 0x00000002
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.CreateToolhelp32Snapshot.argtypes = [wintypes.DWORD, wintypes.DWORD]
    kernel32.CreateToolhelp32Snapshot.restype = wintypes.HANDLE
    kernel32.Process32FirstW.argtypes = [wintypes.HANDLE, ctypes.POINTER(ProcessEntry)]
    kernel32.Process32FirstW.restype = wintypes.BOOL
    kernel32.Process32NextW.argtypes = [wintypes.HANDLE, ctypes.POINTER(ProcessEntry)]
    kernel32.Process32NextW.restype = wintypes.BOOL
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]

    snapshot = kernel32.CreateToolhelp32Snapshot(snapshot_flag, 0)
    invalid_handle = ctypes.c_void_p(-1).value
    if snapshot == invalid_handle:
        return {}
    result: dict[int, tuple[int, str]] = {}
    try:
        entry = ProcessEntry()
        entry.size = ctypes.sizeof(entry)
        success = kernel32.Process32FirstW(snapshot, ctypes.byref(entry))
        while success:
            result[int(entry.process_id)] = (
                int(entry.parent_process_id), entry.executable.lower()
            )
            success = kernel32.Process32NextW(snapshot, ctypes.byref(entry))
    finally:
        kernel32.CloseHandle(snapshot)
    return result


def _windows_rss(process_id: int) -> int:
    import ctypes
    from ctypes import wintypes

    size_type = ctypes.c_size_t

    class ProcessMemoryCounters(ctypes.Structure):
        _fields_ = [
            ("size", wintypes.DWORD),
            ("page_fault_count", wintypes.DWORD),
            ("peak_working_set_size", size_type),
            ("working_set_size", size_type),
            ("quota_peak_paged_pool_usage", size_type),
            ("quota_paged_pool_usage", size_type),
            ("quota_peak_nonpaged_pool_usage", size_type),
            ("quota_nonpaged_pool_usage", size_type),
            ("pagefile_usage", size_type),
            ("peak_pagefile_usage", size_type),
        ]

    query_limited_information = 0x1000
    virtual_memory_read = 0x0010
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    psapi = ctypes.WinDLL("psapi", use_last_error=True)
    kernel32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]
    psapi.GetProcessMemoryInfo.argtypes = [
        wintypes.HANDLE,
        ctypes.POINTER(ProcessMemoryCounters),
        wintypes.DWORD,
    ]
    psapi.GetProcessMemoryInfo.restype = wintypes.BOOL

    handle = kernel32.OpenProcess(
        query_limited_information | virtual_memory_read, False, process_id
    )
    if not handle:
        return 0
    try:
        counters = ProcessMemoryCounters()
        counters.size = ctypes.sizeof(counters)
        if not psapi.GetProcessMemoryInfo(
            handle, ctypes.byref(counters), counters.size
        ):
            return 0
        return int(counters.working_set_size)
    finally:
        kernel32.CloseHandle(handle)


def _linux_process_table() -> dict[int, tuple[int, str]]:
    result: dict[int, tuple[int, str]] = {}
    proc = Path("/proc")
    if not proc.is_dir():
        return result
    for item in proc.iterdir():
        if not item.name.isdigit():
            continue
        try:
            fields = (item / "stat").read_text(encoding="ascii").split()
            result[int(item.name)] = (
                int(fields[3]),
                (item / "comm").read_text(encoding="utf-8").strip().lower(),
            )
        except (OSError, ValueError, IndexError):
            continue
    return result


def _linux_rss(process_id: int) -> int:
    try:
        resident_pages = int(
            (Path("/proc") / str(process_id) / "statm")
            .read_text(encoding="ascii")
            .split()[1]
        )
        return resident_pages * os.sysconf("SC_PAGE_SIZE")
    except (OSError, ValueError, IndexError):
        return 0


def _sample_tree(
    root: int,
    known: set[int],
    executable_names: set[str],
    excluded_processes: set[int],
) -> int | None:
    if os.name == "nt":
        table = _windows_process_table()
        rss = _windows_rss
    elif Path("/proc").is_dir():
        table = _linux_process_table()
        rss = _linux_rss
    else:
        return None

    known.intersection_update(table)
    if root in table:
        known.add(root)
    changed = True
    while changed:
        changed = False
        for process_id, (parent_id, executable) in table.items():
            is_matching_worker = (
                process_id not in excluded_processes
                and executable in executable_names
            )
            if process_id not in known and (
                parent_id in known or is_matching_worker
            ):
                known.add(process_id)
                changed = True
    return sum(rss(process_id) for process_id in known)


def run_measured(
    command: list[str],
    log_path: Path,
    sample_seconds: float = 0.01,
) -> tuple[int, str, float, int | None]:
    """Run a command and return exit, captured text, wall seconds and peak tree RSS."""
    started = time.perf_counter()
    known: set[int] = set()
    peak_rss: int | None = None
    command_paths = [Path(value) for value in command]
    executable_names = {
        value.name.lower() for value in command_paths if value.is_file()
    }
    if not executable_names:
        executable_names.add(command_paths[0].name.lower())
    if os.name == "nt":
        initial_table = _windows_process_table()
    elif Path("/proc").is_dir():
        initial_table = _linux_process_table()
    else:
        initial_table = {}
    excluded_processes = {
        process_id
        for process_id, (_, executable) in initial_table.items()
        if executable in executable_names
    }
    with log_path.open("w", encoding="utf-8") as log:
        process = subprocess.Popen(
            command,
            text=True,
            stdout=log,
            stderr=subprocess.STDOUT,
        )
        while process.poll() is None:
            current = _sample_tree(
                process.pid, known, executable_names, excluded_processes
            )
            if current is not None:
                peak_rss = max(peak_rss or 0, current)
            time.sleep(sample_seconds)
        return_code = process.wait()
        current = _sample_tree(
            process.pid, known, executable_names, excluded_processes
        )
        if current is not None:
            peak_rss = max(peak_rss or 0, current)
    elapsed = time.perf_counter() - started
    return return_code, log_path.read_text(encoding="utf-8"), elapsed, peak_rss
