import sympy
from sympy import symbols, Poly, GF, Matrix, gcd, simplify, cancel
from sympy.polys.domains import GF
from sympy.polys.galoistools import gf_irreducible
import random as pyrandom
from sage.all import *
import time
import matplotlib.pyplot as plt
import networkx as nx
import random


def generate_random_sl2(a_bound: int, b_bound: int, c_bound: int, max_trials: int = 1000):
    """
    Generate a random SL(2, Z) matrix of the form:
        [[1 + 4a, 2b],
         [2c,   1 + d]]

    The parameters a, b, c are chosen uniformly at random from
    [-a_bound, a_bound], [-b_bound, b_bound], [-c_bound, c_bound].
    The integer d is computed so that the determinant = 1.

    Args:
        a_bound (int): maximum absolute value for a.
        b_bound (int): maximum absolute value for b.
        c_bound (int): maximum absolute value for c.
        max_trials (int): maximum number of random attempts.

    Returns:
        tuple: (matrix, (a, b, c, d)) where matrix is a 2x2 list representing the SL2(Z) matrix.

    Raises:
        ValueError: if no valid SL2(Z) matrix is found within max_trials.
    """
    for _ in range(max_trials):
        a = random.randint(-a_bound, a_bound)
        b = random.randint(-b_bound, b_bound)
        c = random.randint(-c_bound, c_bound)

        # Compute d to satisfy det([[1+4a, 2b], [2c, 1+d]]) = 1
        denom = 1 + 4 * a
        if denom == 0:
            continue

        numerator = 4 * (b * c - a)
        if numerator % denom != 0:
            continue

        d = numerator // denom
        mat = [[1 + 4 * a, 2 * b],
               [2 * c,      1 + d]]
        return mat, (a, b, c, d)

    raise ValueError(f"No SL(2,Z) matrix found within {max_trials} trials with given bounds.")


if __name__ == '__main__':
    # Example usage with bounds for a, b, c
    M, params = generate_random_sl2(1000, 1000, 1000)
    print("Generated SL(2,Z) matrix:", M)
    print("Parameters (a, b, c, d):", params)
