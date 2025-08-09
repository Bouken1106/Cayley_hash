import sympy
from sympy import symbols, Poly, GF, Matrix, gcd, simplify, cancel
from sympy.polys.domains import GF
from sympy.polys.galoistools import gf_irreducible
import random as pyrandom
from sage.all import *
import time
import matplotlib.pyplot as plt
import networkx as nx
from sage.all import Matrix, ZZ, random_prime, legendre_symbol, identity_matrix
import random

def random_n_bit_prime(n):
    """
    ビット長 n のランダム素数を返す。
    2^(n-1) ≤ p < 2^n を満たす素数 p を返す。
    """
    lower = 2**(n-1)
    upper = 2**n
    # [lower, upper) の範囲でランダムに素数を返す
    return random_prime(upper, lbound=lower)

# --- 使用例 ---
if __name__ == "__main__":
    for bits in [16, 32, 64, 128]:
        p = random_n_bit_prime(bits)
        print(f"{bits}-bit prime:", p)
        # 実際にビット長を確認
        print(" bit_length =", p.nbits(), "\n")