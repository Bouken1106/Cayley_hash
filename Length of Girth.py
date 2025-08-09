#!/usr/bin/env python3
"""
Compute the girth of the Cayley graph Cay(SL_2(F_p), {A(n),B(n)})
using a simple breadth‑first search (BFS).

Usage (example):
    $ python girth_bfs.py 101 2 25
→ explores up to depth 25 (cycle length 50) for p=101, n=2.

If SageMath is available, uncomment the [SAGE] section to accelerate
matrix operations.  Both pure‑Python and Sage paths share the same API.
"""

import sys
import math
from collections import deque
from typing import Tuple, Dict, Optional
import RandomPrime

# -------------------- matrix helpers --------------------
# Pure‑Python representation: 4‑tuple (a,b,c,d) in row‑major order.
Matrix = Tuple[int, int, int, int]

def mat_mul(x: Matrix, y: Matrix, p: int) -> Matrix:
    a, b, c, d = x
    e, f, g, h = y
    return (
        (a * e + b * g) % p,
        (a * f + b * h) % p,
        (c * e + d * g) % p,
        (c * f + d * h) % p,
    )

def mat_inv(m: Matrix, p: int) -> Matrix:
    a, b, c, d = m  # det == 1 なので逆行列は (d,-b,-c,a)
    return (d % p, (-b) % p, (-c) % p, a % p)

# -------------------- BFS core --------------------

def girth_bfs(p: int, n: int, max_depth: int = 25) -> Optional[int]:
    """Return the girth (even integer) if a cycle ≤ 2*max_depth exists; else None."""

    # generators A(n), B(n) and their inverses
    A: Matrix = (1 % p, n % p, 0, 1 % p)
    B: Matrix = (1 % p, 0, n % p, 1 % p)
    gens = (A, mat_inv(A, p), B, mat_inv(B, p))

    I: Matrix = (1, 0, 0, 1)

    dist: Dict[Matrix, int] = {I: 0}
    parent: Dict[Matrix, Matrix] = {I: I}
    q = deque([I])
    min_cycle = math.inf

    while q:
        g = q.popleft()
        d = dist[g]
        if d * 2 + 1 >= min_cycle:
            continue  # cannot find shorter cycles beyond this depth
        if d >= max_depth:
            continue  # depth limit reached

        for h in gens:
            nxt = mat_mul(h, g, p)
            if nxt not in dist:
                dist[nxt] = d + 1
                parent[nxt] = g
                q.append(nxt)
            elif parent[g] != nxt:  # non‑tree edge → cycle detected
                cycle_len = d + dist[nxt] + 1
                if cycle_len < min_cycle:
                    min_cycle = cycle_len

    return None if min_cycle is math.inf else min_cycle

# -------------------- CLI --------------------

def _parse_int(s: str, name: str) -> int:
    try:
        v = int(s)
        if v <= 0:
            raise ValueError
        return v
    except ValueError:
        sys.exit(f"{name} must be a positive integer: {s}")

if __name__ == "__main__":
    if len(sys.argv) not in {1, 4}:
        print(
            "Usage: python girth_bfs.py [p n max_depth]\n"
            "  p          : prime modulus (default 101)\n"
            "  n          : parameter for A(n),B(n)  (default 2)\n"
            "  max_depth  : search up to cycle length 2*max_depth (default 25)",
            file=sys.stderr,
        )
        sys.exit(1)

    p = _parse_int(sys.argv[1], "p") if len(sys.argv) >= 2 else random_n_bit_prime(10)
    n = _parse_int(sys.argv[2], "n") if len(sys.argv) >= 3 else 2
    depth = _parse_int(sys.argv[3], "max_depth") if len(sys.argv) == 4 else 25

    g = girth_bfs(p, n, depth)
    if g is None:
        print(f"No cycle found up to length {2 * depth}.")
    else:
        print(f"girth = {g}")
