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





def lift_identity(p, c_bound=20, trials=4000):
    """
    Tillich–Zémor『Group-theoretic hash functions』Algorithm 3 に基づき，
    素数 p について次を満たす行列 M ∈ SL₂(ℤ) を返す：
        • M ≡ I (mod p)
        • det M = 1
        • M ≠ I   （非自明リフト）

    Parameters
    ----------
    p : int  または Sage Integer
        既知の素数
    c_bound : int
        試行する c を 1 … c_bound からランダムに選ぶ
    trials : int
        (c, p′) の組合せを試す回数
    """
                           
    bitlen = p.bit_length()         # Python/Sage 共通のビット長取得

    for _ in range(trials):
        # 1) 小さな整数 c をランダムに選択
        c = random.randint(1, c_bound)

        # 2) p と同程度のビット長を持つランダム素数 p′ を選択
        p_prime = random_prime(2**bitlen, lbound=max(2**(bitlen - 1),3))

        # 3) Δ = c² p² + 4c が p′ 上で平方剰余か判定
        Delta = c**2 * p**2 + 4 * c
        if legendre_symbol(Delta, p_prime) != 1:
            continue  # 非平方 → 次の試行

        # 4) Δ の平方根 (mod p′) をすべて取得
        roots = list(Mod(Delta, p_prime).sqrt(all=True))
        if not roots:
            continue

        for root_mod in roots:          # root_mod は Z/p′Z 上の元
            root = int(root_mod)        # Python int へキャスト
            # k₁ = (c·p + root)/2 が整数か？
            if (c * p + root) % 2 != 0:
                continue
            k1 = (c * p + root) // 2
            k4 = c * p - k1

            # k₂ = (c + c·p·k₁ − k₁²) / p′ が整数か？
            k2_num = c + c * p * k1 - k1**2
            if k2_num % p_prime != 0:
                continue
            k2 = k2_num // p_prime
            k3 = p_prime

            # 行列 M を構成
            M = Matrix(ZZ, [
                [1 + k1 * p, k2 * p],
                [k3 * p,     1 + k4 * p]
            ])

            # 条件検証
            if M.det() == 1 and M != identity_matrix(ZZ, 2):
                return M

    raise ValueError(f"非自明リフトが見つかりませんでした (p = {p})")

# --- 使用例 ---
if __name__ == "__main__":
    p = 3
    M = lift_identity(p)
    print("Lifted matrix M =\n", M)
    print("\nM mod p =\n", M.apply_map(lambda x: x % p))
    print("\nM mod p =\n", M.apply_map(lambda x: x % 4))
