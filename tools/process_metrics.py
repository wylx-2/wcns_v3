#!/usr/bin/env python3
"""Small dependency-free process-tree RSS sampler for release measurements."""

from __future__ import annotations

import os
import subprocess
import threading
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
            result[int(entry.process_id)] = (int(entry.parent_process_id), entry.executable.lower())
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
        if not psapi.GetProcessMemoryInfo(handle, ctypes.byref(counters), counters.size):
            return 0
        return int(counters.working_set_size)
    finally:
        kernel32.CloseHandle(handle)


def _windows_set_high_priority(process_id: int) -> bool:
    import ctypes
    from ctypes import wintypes

    process_set_information = 0x0200
    query_limited_information = 0x1000
    high_priority_class = 0x00000080
    kernel32 = ctypes.WinDLL("kernel32", use_last_error=True)
    kernel32.OpenProcess.argtypes = [wintypes.DWORD, wintypes.BOOL, wintypes.DWORD]
    kernel32.OpenProcess.restype = wintypes.HANDLE
    kernel32.SetPriorityClass.argtypes = [wintypes.HANDLE, wintypes.DWORD]
    kernel32.SetPriorityClass.restype = wintypes.BOOL
    kernel32.CloseHandle.argtypes = [wintypes.HANDLE]

    handle = kernel32.OpenProcess(
        process_set_information | query_limited_information, False, process_id
    )
    if not handle:
        return False
    try:
        return bool(kernel32.SetPriorityClass(handle, high_priority_class))
    finally:
        kernel32.CloseHandle(handle)


def _raise_new_windows_workers(
    executable_names: set[str],
    excluded_processes: set[int],
    expected_processes: int,
) -> None:
    deadline = time.perf_counter() + 2.0
    raised: set[int] = set()
    while time.perf_counter() < deadline:
        table = _windows_process_table()
        candidates = {
            process_id
            for process_id, (_, executable) in table.items()
            if executable in executable_names and process_id not in excluded_processes
        }
        for process_id in candidates - raised:
            if not _windows_set_high_priority(process_id):
                raise RuntimeError(f"could not set high priority for process {process_id}")
            raised.add(process_id)
        if len(raised) >= expected_processes:
            return
        time.sleep(0.01)
    raise RuntimeError(
        "could not find all performance workers before the priority deadline: "
        f"expected {expected_processes}, found {len(raised)}"
    )


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
            (Path("/proc") / str(process_id) / "statm").read_text(encoding="ascii").split()[1]
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
                process_id not in excluded_processes and executable in executable_names
            )
            if process_id not in known and (parent_id in known or is_matching_worker):
                known.add(process_id)
                changed = True
    return sum(rss(process_id) for process_id in known)


def run_measured(
    command: list[str],
    log_path: Path,
    sample_seconds: float | None = 0.01,
    high_priority_names: set[str] | None = None,
    expected_high_priority_processes: int = 0,
) -> tuple[int, str, float, int | None]:
    """Run a command and return exit, captured text, wall seconds and peak tree RSS."""
    started = time.perf_counter()
    known: set[int] = set()
    peak_rss: int | None = None
    command_paths = [Path(value) for value in command]
    executable_names = {value.name.lower() for value in command_paths if value.is_file()}
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
        if high_priority_names:
            if os.name != "nt":
                raise RuntimeError("high-priority worker selection is Windows-only")
            if expected_high_priority_processes <= 0:
                process.terminate()
                process.wait()
                raise RuntimeError("expected high-priority process count must be positive")
            try:
                _raise_new_windows_workers(
                    {value.lower() for value in high_priority_names},
                    excluded_processes,
                    expected_high_priority_processes,
                )
            except BaseException:
                process.terminate()
                process.wait()
                raise
        sampler_error: list[BaseException] = []
        stop_sampler = threading.Event()
        sampler: threading.Thread | None = None
        if sample_seconds is not None:

            def sample_until_stopped() -> None:
                nonlocal peak_rss
                try:
                    while not stop_sampler.is_set():
                        current = _sample_tree(
                            process.pid, known, executable_names, excluded_processes
                        )
                        if current is not None:
                            peak_rss = max(peak_rss or 0, current)
                        stop_sampler.wait(sample_seconds)
                except BaseException as error:
                    sampler_error.append(error)

            sampler = threading.Thread(target=sample_until_stopped, daemon=True)
            sampler.start()
        return_code = process.wait()
        finished = time.perf_counter()
        stop_sampler.set()
        if sampler is not None:
            sampler.join()
        if sampler_error:
            raise RuntimeError("process sampler failed") from sampler_error[0]
        if sample_seconds is not None:
            current = _sample_tree(process.pid, known, executable_names, excluded_processes)
            if current is not None:
                peak_rss = max(peak_rss or 0, current)
    elapsed = finished - started
    return return_code, log_path.read_text(encoding="utf-8"), elapsed, peak_rss
