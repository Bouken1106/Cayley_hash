#!/usr/bin/env python3
"""
Factorization in the subgroup H = <A(2), B(2)> of SL_2(Z).

A(2) = [[1, 2], [0, 1]]
B(2) = [[1, 0], [2, 1]]

Given a matrix M ∈ H, the function ``factorize`` returns a list of generators
['A', 'A^-1', 'B', 'B^-1'] whose product equals M (multiplication is read *left
→ right* exactly as you would write a word in a free group).

Implementation follows the reduction algorithm provided by the user.  One subtle
point: after we have reduced (left‑multiplied) M to the identity using a word
S = g_k … g_1, we obtain S·M = I and hence M = g_1⁻¹ … g_k⁻¹.  **Note that the
order is *not* reversed** – we only invert every generator.  (My first draft had
an unnecessary ``reversed``; this version fixes the bug.)
"""
from __future__ import annotations
from typing import List
import random
# ---------------------------------------------------------------------------
# low‑level helpers
Matrix = List[List[int]]

def mul(X: Matrix, Y: Matrix) -> Matrix:
    """Return X·Y for 2×2 integer matrices (arbitrary precision)."""
    return [
        [X[0][0] * Y[0][0] + X[0][1] * Y[1][0],
         X[0][0] * Y[0][1] + X[0][1] * Y[1][1]],
        [X[1][0] * Y[0][0] + X[1][1] * Y[1][0],
         X[1][0] * Y[0][1] + X[1][1] * Y[1][1]]
    ]

def det(M: Matrix) -> int:
    return M[0][0] * M[1][1] - M[0][1] * M[1][0]

# ---------------------------------------------------------------------------
# generators
A      : Matrix = [[1, 2], [0, 1]]
A_inv  : Matrix = [[1, -2], [0, 1]]
B      : Matrix = [[1, 0], [2, 1]]
B_inv  : Matrix = [[1, 0], [-2, 1]]
GEN = {'A': A, 'A^-1': A_inv, 'B': B, 'B^-1': B_inv}
INV = {'A': 'A^-1', 'A^-1': 'A', 'B': 'B^-1', 'B^-1': 'B'}

# ---------------------------------------------------------------------------
# main algorithm

def factorize(M: Matrix, max_steps: int = 10_000) -> List[str]:
    """Return a word in {A(2)±¹, B(2)±¹} whose product equals *M*.

    The algorithm is deterministic and terminates for every element of H,
    strictly decreasing |b|+|c| at each step.
    """
    if det(M) != 1:
        raise ValueError("input matrix must have determinant 1")

    W = [row[:] for row in M]        # working copy
    ops: List[str] = []             # record left multiplications

    for _ in range(max_steps):
        if W == [[1, 0], [0, 1]]:
            break
        a, b = W[0]
        c, d = W[1]

        if abs(b) > abs(d):            # reduce |b|
            gen = 'A^-1' if b * d > 0 else 'A'
        elif abs(c) > abs(a):          # reduce |c|
            gen = 'B^-1' if a * c > 0 else 'B'
        else:                          # should not happen
            raise RuntimeError("algorithm stalled — is the matrix in H?")

        W = mul(GEN[gen], W)           # left multiply
        ops.append(gen)
    else:
        raise RuntimeError("exceeded max_steps, possible invalid input")

    # We have (g_k … g_1)·M = I ⇒ M = g_1⁻¹ … g_k⁻¹ (no reversal!)
    return [INV[g] for g in ops]

def random_word(n):
    """
    長さ n のランダムな word を生成し、リストで返します。
    要素は 'A', 'B', 'A^-1', 'B^-1' のいずれかです。
    """
    generators = ['A', 'B', 'A^-1', 'B^-1']
    return [random.choice(generators) for _ in range(n)]


# ---------------------------------------------------------------------------
# demo / self‑test
if __name__ == "__main__":
    from functools import reduce
    word_size = 100
    word = random_word(word_size)
    print('Input word:', word)

    # build M = A^3 · B · A^-2 for a quick test
    example = word
    M = reduce(mul, (GEN[g] for g in example))
    print("Input matrix :", M)

    word = factorize(M)
    print("Recovered    :", " ".join(word))

    M_rec = reduce(mul, (GEN[g] for g in word))
    assert M_rec == M, "factorization failed!"
    print("✔ Verification passed")

    word = factorize([[-7, 956], [110, -15023]])
    print("Recovered    :", " ".join(word))
