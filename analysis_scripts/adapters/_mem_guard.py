"""Shared 55 GiB process budget for Windows benchmark adapters.

Two layers enforce the same limit for Python and wrapped R processes:

1. require_under(est_gb): pre-flight. Estimate peak RAM from data dimensions and
   sys.exit(1) BEFORE allocating if it exceeds the budget. Nothing gets allocated.
2. cap_self(gb): hard backstop. Put the current process (and the worker children it
   spawns, which inherit the job) into a Windows Job Object with a committed-memory
   limit. Exceeding the limit terminates the child with a non-zero status.

CLI launcher mode (for wrapping a child exe, e.g. Rscript):
    python _mem_guard.py --limit-gb 55 -- <exe> <args...>
caps itself then runs <exe> as a child (which inherits the capped job) and forwards
the child's exit code.
"""
import ctypes
import sys
from ctypes import wintypes

DEFAULT_BUDGET_GB = 55

# --- Windows Job Object plumbing (pure ctypes, no extra dependency) ----------
_JobObjectExtendedLimitInformation = 9
_JOB_OBJECT_LIMIT_JOB_MEMORY = 0x00000200
_ULONG_PTR = ctypes.c_size_t
_HELD = []  # keep job handles alive for the process lifetime


class _IO_COUNTERS(ctypes.Structure):
    _fields_ = [(n, ctypes.c_ulonglong) for n in (
        "ReadOperationCount", "WriteOperationCount", "OtherOperationCount",
        "ReadTransferCount", "WriteTransferCount", "OtherTransferCount")]


class _BASIC_LIMIT(ctypes.Structure):
    _fields_ = [
        ("PerProcessUserTimeLimit", ctypes.c_longlong),
        ("PerJobUserTimeLimit", ctypes.c_longlong),
        ("LimitFlags", wintypes.DWORD),
        ("MinimumWorkingSetSize", _ULONG_PTR),
        ("MaximumWorkingSetSize", _ULONG_PTR),
        ("ActiveProcessLimit", wintypes.DWORD),
        ("Affinity", _ULONG_PTR),
        ("PriorityClass", wintypes.DWORD),
        ("SchedulingClass", wintypes.DWORD),
    ]


class _EXTENDED_LIMIT(ctypes.Structure):
    _fields_ = [
        ("BasicLimitInformation", _BASIC_LIMIT),
        ("IoInfo", _IO_COUNTERS),
        ("ProcessMemoryLimit", _ULONG_PTR),
        ("JobMemoryLimit", _ULONG_PTR),
        ("PeakProcessMemoryUsed", _ULONG_PTR),
        ("PeakJobMemoryUsed", _ULONG_PTR),
    ]


def cap_self(gb=DEFAULT_BUDGET_GB):
    """Best-effort hard cap on committed memory for this process + its children."""
    if sys.platform != "win32":
        return False
    k = ctypes.WinDLL("kernel32", use_last_error=True)
    k.CreateJobObjectW.restype = wintypes.HANDLE
    k.CreateJobObjectW.argtypes = [wintypes.LPVOID, wintypes.LPCWSTR]
    k.SetInformationJobObject.restype = wintypes.BOOL
    k.SetInformationJobObject.argtypes = [wintypes.HANDLE, ctypes.c_int, wintypes.LPVOID, wintypes.DWORD]
    k.GetCurrentProcess.restype = wintypes.HANDLE
    k.AssignProcessToJobObject.restype = wintypes.BOOL
    k.AssignProcessToJobObject.argtypes = [wintypes.HANDLE, wintypes.HANDLE]

    h = k.CreateJobObjectW(None, None)
    if not h:
        return False
    info = _EXTENDED_LIMIT()
    info.BasicLimitInformation.LimitFlags = _JOB_OBJECT_LIMIT_JOB_MEMORY
    info.JobMemoryLimit = int(gb * (1024 ** 3))
    if not k.SetInformationJobObject(h, _JobObjectExtendedLimitInformation,
                                     ctypes.byref(info), ctypes.sizeof(info)):
        return False
    if not k.AssignProcessToJobObject(h, k.GetCurrentProcess()):
        return False
    _HELD.append(h)
    return True


def require_under(est_gb, limit_gb=DEFAULT_BUDGET_GB, what=""):
    """Pre-flight: fail before allocating if the estimated peak exceeds the budget."""
    if est_gb is not None and est_gb > limit_gb:
        sys.stderr.write(
            "FAIL: estimated peak memory %.1f GB exceeds %g GB budget for %s.\n"
            % (est_gb, limit_gb, what))
        sys.exit(1)


def dense_gb(n_cells, n_genes, bytes_per=8, copies=1):
    return float(n_cells) * float(n_genes) * bytes_per * copies / (1024 ** 3)


def square_gb(n, bytes_per=8, copies=1):
    return float(n) * float(n) * bytes_per * copies / (1024 ** 3)


def _cli(argv):
    if "--" not in argv:
        sys.stderr.write("usage: _mem_guard.py --limit-gb N -- <exe> <args...>\n")
        return 2
    limit = DEFAULT_BUDGET_GB
    if "--limit-gb" in argv:
        limit = float(argv[argv.index("--limit-gb") + 1])
    cmd = argv[argv.index("--") + 1:]
    if not cmd:
        sys.stderr.write("_mem_guard.py: no command after --\n")
        return 2
    cap_self(limit)  # child inherits this job -> child is capped too
    import subprocess
    return subprocess.run(cmd).returncode


if __name__ == "__main__":
    sys.exit(_cli(sys.argv[1:]))
