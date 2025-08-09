import sympy
from sympy import symbols, Poly, GF, Matrix, gcd, simplify, cancel
from sympy.polys.domains import GF
from sympy.polys.galoistools import gf_irreducible
import random as pyrandom
from sage.all import *
from sage.all import GF, Matrix
import time
import matplotlib.pyplot as plt
import networkx as nx
from collections import deque

from typing import List, Tuple, Union

# Define the standard generators A(2), B(2) for SL2(Z)
A2 = [[1, 2],
      [0, 1]]
B2 = [[1, 0],
      [2, 1]]

# Their inverses in SL2(Z)
# For A2 = [[1,2],[0,1]], A2^{-1} = [[1,-2],[0,1]]
# For B2 = [[1,0],[2,1]], B2^{-1} = [[1,0],[-2,1]]
A2_inv = [[1, -2],
          [0,  1]]
B2_inv = [[1,  0],
          [-2, 1]]

# Helper: multiply two 2x2 matrices
def mat_mul(X: List[List[int]], Y: List[List[int]]) -> List[List[int]]:
    return [
        [X[0][0]*Y[0][0] + X[0][1]*Y[1][0],  X[0][0]*Y[0][1] + X[0][1]*Y[1][1]],
        [X[1][0]*Y[0][0] + X[1][1]*Y[1][0],  X[1][0]*Y[0][1] + X[1][1]*Y[1][1]],
    ]

# Compute the sum of absolute values of entries of a 2x2 matrix
def abs_sum(M: List[List[int]]) -> int:
    return sum(abs(M[i][j]) for i in range(2) for j in range(2))

# Check if a matrix is the identity
def is_identity(M: List[List[int]]) -> bool:
    return M == [[1,0],[0,1]]

# Decompose M into a word in A2, B2 and their inverses by greedy reduction
# Returns either a list of operations or False
# Each operation is a tuple (side, gen_name), where side is 'L' or 'R'
# and gen_name is one of 'A2', 'B2', 'A2_inv', 'B2_inv'.
def decompose_sl2(M: List[List[int]]) -> Union[List[Tuple[str,str]], bool]:
    # Current matrix copy
    current = [row[:] for row in M]
    # Record of operations
    ops: List[Tuple[str,str]] = []
    
    while True:
        curr_sum = abs_sum(current)
        best_sum = curr_sum
        best_op = None
        best_mat = None
        # Try 8 possible moves
        for side in ('L', 'R'):
            for gen_name, gen_mat in (
                ('A2', A2), ('B2', B2), ('A2_inv', A2_inv), ('B2_inv', B2_inv)
            ):
                if side == 'L':
                    candidate = mat_mul(gen_mat, current)
                else:
                    candidate = mat_mul(current, gen_mat)
                s = abs_sum(candidate)
                if s < best_sum:
                    best_sum = s
                    best_op = (side, gen_name)
                    best_mat = candidate
        # If no reduction found, stop
        if best_op is None:
            break
        # Apply the best operation
        current = best_mat  # type: ignore
        ops.append(best_op)
    
    # Check if reduced to identity
    if is_identity(current):
        return ops
    else:
        return False

# Example usage:
if __name__ == '__main__':
    # Example matrix in SL2(Z)
    M = [[5, 2], [2, 1]]  # this is C(2) = A2*B2
    result = decompose_sl2(M)
    if result is False:
        print("No decomposition found.")
    else:
        print("Decomposition operations (in order applied):")
        print(result)
