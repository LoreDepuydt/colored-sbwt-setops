#!/usr/bin/env python3
# Prints how many bytes of the given files are in the page cache, and their total size (as fincore,
# which not every server has): "<cached_bytes> <total_bytes>". Uses mincore(2) on a read-only mapping.
# Usage: python3 -I cached_bytes.py <file>...
import ctypes, mmap, os, sys

libc = ctypes.CDLL(None, use_errno=True)
libc.mmap.restype = ctypes.c_void_p
libc.mmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_long]
libc.munmap.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
libc.mincore.argtypes = [ctypes.c_void_p, ctypes.c_size_t, ctypes.c_char_p]
MAP_FAILED = ctypes.c_void_p(-1).value
page = mmap.PAGESIZE
cached = total = 0
for path in sys.argv[1:]:
    size = os.path.getsize(path)
    total += size
    if size == 0:
        continue
    with open(path, 'rb') as f:
        a = libc.mmap(None, size, mmap.PROT_READ, mmap.MAP_SHARED, f.fileno(), 0)
        if a in (None, MAP_FAILED):
            sys.exit(f"mmap failed on {path}: {os.strerror(ctypes.get_errno())}")
        vec = ctypes.create_string_buffer((size + page - 1) // page)
        ok = libc.mincore(a, size, vec) == 0
        err = ctypes.get_errno()
        libc.munmap(a, size)
        if not ok:
            sys.exit(f"mincore failed on {path}: {os.strerror(err)}")
        cached += min(sum(b & 1 for b in vec.raw) * page, size)
print(cached, total)
