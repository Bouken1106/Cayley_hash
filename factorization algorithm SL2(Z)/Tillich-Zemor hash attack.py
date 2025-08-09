import sympy
from sympy import symbols, Poly, GF, Matrix, gcd
from sympy.polys.domains import GF
from sympy.polys.galoistools import gf_irreducible
import random as pyrandom
from sage.all import *
import time
import matplotlib.pyplot as plt


# 変数 x の定義
x = symbols('x')


# GF(2)上の次数3既約多項式を1つランダムに生成
def random_irreducible_polynomial(degree):
    F = GF(2)                          
    R = PolynomialRing(F, 'x')         
    x_poly = R.gen()                     
    while True:
        f_candidate = R.random_element(degree)
        if f_candidate.is_irreducible():
            f = f_candidate
            break
    return f


def compute_p_n_minus_1(q_poly):
    """
    GF(2) 上で、行列 A * u = b の一意解（自由変数は 0 とする）を
    ガウス・ジョルダン法により求める。
    A は m×n Matrix、b は長さ m のリスト（各要素は 0 または 1）。
    """
    def solve_mod2(A, b):
        m, n = A.nrows(), A.ncols()
        # A と b を拡大係数行列 M にまとめる（各成分を整数 0,1 とする）
        M = []
        for i in range(m):
            row = [int(A[i, j] % 2) for j in range(n)]
            row.append(int(b[i] % 2))
            M.append(row)
            
        pivot_row = 0
        pivot_cols = []
        for j in range(n):
            # 列 j で 1 のピボットを探す
            pivot = None
            for i in range(pivot_row, m):
                if M[i][j] == 1:
                    pivot = i
                    break
            if pivot is None:
                continue
            # ピボット行と pivot_row の行を入れ替え
            if pivot != pivot_row:
                M[pivot_row], M[pivot] = M[pivot], M[pivot_row]
            pivot_cols.append(j)
            # ピボット行以外のすべての行について、列 j を消去
            for i in range(m):
                if i != pivot_row and M[i][j] == 1:
                    # 行の加算（GF(2) では XOR 相当）
                    M[i] = [(M[i][k] + M[pivot_row][k]) % 2 for k in range(n+1)]
            pivot_row += 1
            if pivot_row == m:
                break

        # 一致しない方程式のチェック
        for i in range(m):
            if all(M[i][j] == 0 for j in range(n)) and M[i][n] != 0:
                raise ValueError("No solution exists in GF(2).")
        # 解 u はピボット列に沿って（自由変数は 0 とする）決定される。
        u = [0] * n
        for i, col in enumerate(pivot_cols):
            u[col] = M[i][n]
        return u
    """
    入力: q_poly -- GF(2) 上の既約多項式（sympy Poly, modulus=2）
    出力: p_poly -- Mesirov–Sweet の手法で得られる p_(n-1)(x)（Poly, modulus=2）
           u      -- 線形系の解 (リスト)
           A      -- 係数行列 (sympy Matrix)
           b      -- 右辺ベクトル (list)
    """
    n = q_poly.degree()
    
    # polynomials g_i, i=0,...,n を構築。各 g_i を q_poly で剰余計算して次数 < n にする。
    g_list = []
    # g_0 = 1
    g0 = Poly(1, x, modulus=2)  # 修正: modulus=2 を使用
    g_list.append(g0)
    
    for i in range(1, n+1):
        # g_i = x^(i-1) + x^(2*i-1) + x^(2*i)
        expr = x**(i-1) + x**(2*i - 1) + x**(2*i)
        gi = Poly(expr, x, modulus=2).rem(q_poly)
        # 確実に次数 < n となるよう変換
        gi = Poly(gi.as_expr(), x, modulus=2)  # 修正: modulus=2 を使用
        g_list.append(gi)
    
    # 行列 A: 各 g_i について、係数 (x^0～x^(n-1)) を抽出し、(n+1)×n 行列を構築
    A_rows = []
    for gi in g_list:
        # 係数を 0～n-1 次でリスト化
        coeffs = [int(gi.coeff_monomial(x**j)) % 2 for j in range(n)]
        A_rows.append(coeffs)
    A = Matrix(A_rows)
    
    # 右辺ベクトル b: 長さ n+1、最初と最後が 1、その他は 0
    b = [1] + [0]*(n-1) + [1]
    
    # GF(2) 上で A * u = b を解く（自由変数は 0 とする）
    print("solving equation")
    u = solve_mod2(A, b)
    
    # Laurent 多項式 L(x)= sum_{i=1}^n u_i * x^(-i) を考え、
    # p(x) = 非負冪項（x^0以上）の部分＝ sum_{i=1}^{n} sum_{j=i}^{n} u_i * (qの x^j の係数) * x^(j-i)
    q_coeffs = [int(q_poly.coeff_monomial(x**j)) % 2 for j in range(n+1)]  # j=0..n
    # p(x) の係数（次数は 0～n-1 となるはず）
    p_coeffs = [0] * n
    for i in tqdm(range(1, n+1)):
        if u[i-1] == 0:
            continue
        for j in range(i, n+1):  # j = i, i+1, …, n
            exp = j - i  # 対応する x の指数
            p_coeffs[exp] = (p_coeffs[exp] + q_coeffs[j]) % 2

    # p(x) の構築
    p_expr = sum(p_coeffs[k]*x**k for k in range(n))
    p_poly = Poly(p_expr, x, modulus=2)  # 修正: modulus=2 を使用
    
    return p_poly

#qから得られたp_(n-1)に対して衝突列を構成
def get_collision_bit_sequence(q_poly, p_n_minus_1_poly):
    """
    q_poly: 既約多項式 q(x) = p_n(x) (degree n, GF(2))
    p_n_minus_1_poly: p_(n-1)(x) (degree n-1, GF(2))
    
    q_poly と p_n_minus_1_poly から、Euclidean 互除法を実行して、
    各ステップの商が形 x + β_i (β_i ∈ {0,1}) となると仮定し、
    その商の定数項を順に収集します。
    
    ※ ここで得られた商の順序は
        [β_n, β_(n-1), ..., β_1] となるので、反転して
        [β_1, β_2, ..., β_n] を得ます。
    ※ 論文の通り、p1 = x + b1 + 1 から b1 = β_1+1 (mod 2) となり、
        i > 1 の場合は b_i = β_i とします。
    
    その後、得られたビット列 (b_1, ..., b_n) を元に、衝突に用いる回文
    v = b_n b_(n-1) ... b_1 b_1 ... b_(n-1) b_n
    （文字列として連結）を生成して返します。
    
    戻り値: 
      - bit_list: [b_1, b_2, ..., b_n]（リスト形式）
      - collision_bit_string: 回文ビット列の文字列（例: "1011101110"）
    """
    # Euclidean 互除法による連分数展開
    a = q_poly
    b = p_n_minus_1_poly
    beta_list = []  # 各商の定数項を格納（順序は最初の商が β_n となる）
    while b.degree() >= 0 and b != 0:
        quo, rem = a.div(b)
        if quo.degree() != 1:
            raise ValueError("商の次数が 1 であることを期待しましたが、得られた次数は {}".format(quo.degree()))
        # quo は形 x + c として、c を取得（c = β_i）
        c = int(quo.coeff_monomial(x**0)) % 2
        beta_list.append(c)
        a = b
        b = rem
    # beta_list の順序は [β_n, β_(n-1), ..., β_1] であるため、反転する
    beta_list_reversed = beta_list[::-1]  # [β_1, β_2, ..., β_n]
    
    # 論文の定義: b1 = β_1+1 (mod 2), それ以外 b_i = β_i (i>1)
    b_list = []
    for i, beta in enumerate(beta_list_reversed):
        if i == 0:
            b_list.append((beta + 1) % 2)
        else:
            b_list.append(beta)
    
    # 衝突に用いる回文の生成:
    # b_list = [b_1, b_2, ..., b_n] から
    # v = b_n ... b_1 b_1 ... b_n を生成
    first_half = b_list[::-1]  # [b_n, ..., b_1]
    second_half = b_list  # [b_1, ..., b_n]
    collision_bits = first_half + second_half  # 回文: 全長 2n
    collision_bit_string = "".join(str(bit) for bit in collision_bits)
    
    return b_list, collision_bit_string
time_list = []
k = 3
l = 3
for i in range(3,l):
    y = 0
    for _ in range(k):
        q_poly = Poly(random_irreducible_polynomial(i),x,modulus=2)
        q_poly.as_expr()

        start = time.perf_counter()

        p_n_minus_1_poly = compute_p_n_minus_1(q_poly)
        p_n_minus_1_poly.as_expr()
        bit_list, collision_str = get_collision_bit_sequence(q_poly, p_n_minus_1_poly)
        
        end = time.perf_counter()
        y = y + end - start
    
    time_list.append(y/k)

plt.plot(list(range(l - 3)), time_list, marker='o')
plt.xlabel("")
plt.ylabel("")
plt.grid(True)
plt.show()

if 1 == 1:
    q_poly = Poly(x**1021 + x**5 + x**2 + x + 1,x,modulus=2)
    q_poly.as_expr()

    p_n_minus_1_poly = compute_p_n_minus_1(q_poly)
    p_n_minus_1_poly.as_expr()
    bit_list, collision_str = get_collision_bit_sequence(q_poly, p_n_minus_1_poly)
    print(bit_list, collision_str)

