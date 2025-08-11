# -*- coding: utf-8 -*-
"""
Sage(Python) 実行向けの BFS による girth 探索コード。

【重要】このファイルは通常の Python ではなく、Sage の Python で実行してください。
    実行例:  sage -python Cayley_girth_BFS_SL2Fp_sage_python.py

通常の Python で動かしたい場合は、GF/matrix を使わない純 Python 版を使ってください（別キャンバスの .py 版）。
"""

# Sage のオブジェクトを Python スクリプトで使うための import
from sage.all import GF, matrix, identity_matrix
from collections import deque
from typing import Optional, List, Tuple

# 型エイリアス（メモ用。Sage の行列型自体は厳密にはこの別名ではない）
Mat = type(matrix(GF(3), [[1,0],[0,1]]))

# ------------------------------
# 生成子 A(n), B(n) とその逆元（GF(p) 上）
# ------------------------------

def generators_A_B(n: int, p: int):
    """GF(p) 上の A(n), B(n) とその逆元を返す。

    A(n)   = [[1, n],[0, 1]]
    B(n)   = [[1, 0],[n, 1]]
    A^{-1} = [[1,-n],[0, 1]]
    B^{-1} = [[1, 0],[-n,1]]

    返り値:
      gens    : (A, B, Ainv, Binv)
      inv_idx : 逆元のインデックス対応（例: inv_idx[0] == 2）
    """
    F = GF(p)
    nF = F(n)
    A    = matrix(F, [[1,  nF],[0, 1]])
    B    = matrix(F, [[1,   0],[nF,1]])
    Ainv = matrix(F, [[1, -nF],[0, 1]])
    Binv = matrix(F, [[1,   0],[-nF,1]])
    gens = (A, B, Ainv, Binv)
    inv_idx = {0:2, 1:3, 2:0, 3:1}
    return gens, inv_idx

# ------------------------------
# キー化（辞書のキーとして高速化: 4成分の整数タプル化）
# ------------------------------

def key_of(M: Mat) -> Tuple[int,int,int,int]:
    """行列 M を辞書キーとして使いやすい 4 つ組（int）に変換。"""
    return (int(M[0,0]), int(M[0,1]), int(M[1,0]), int(M[1,1]))

# ------------------------------
# 語の適用＆表示（検算・デバッグ用）
# ------------------------------

def apply_word(word: List[int], gens) -> Mat:
    """生成子インデックス列 word を左から順に掛けて得られる行列を返す。"""
    F = gens[0].base_ring()
    M = identity_matrix(F, 2)
    for i in word:
        M = M * gens[i]
    return M


def pretty_word(word: List[int]) -> str:
    """インデックス列を可読な記号列に整形（A, B, a, b）。"""
    sym = ['A', 'B', 'a', 'b']
    return ''.join(sym[i] for i in word)

# ------------------------------
# BFS による girth 探索本体
# ------------------------------

def cayley_girth_bfs(p: int, n: int, max_depth: int = 40, return_word: bool = True):
    """Cay(SL_2(F_p), {A(n)^±1, B(n)^±1}) の **girth** を BFS で探索。

    引数:
      p           : 素数（法）
      n           : A(n), B(n) のパラメータ
      max_depth   : 探索の最大深さ（打切り）。深いほど指数的に重い。
      return_word : True なら (girth, word) を返す。False なら girth のみ。

    返り値:
      - return_word=True  : (girth, word_indices) もしくは None
      - return_word=False : girth（int）もしくは None
    """
    F = GF(p)
    I = identity_matrix(F, 2)
    Ikey = key_of(I)

    gens, inv_idx = generators_A_B(n, p)

    # dist: 各状態（行列）への最短到達深さ
    dist = {Ikey: 0}

    # parent: パス復元用（その状態に来る直前のキーと使った生成子 idx）
    parent = {Ikey: (None, None)}

    # キュー: (行列, 最後に使った生成子 idx)。開始時は -1（なし）
    q = deque([(I, -1)])

    while q:
        M, last = q.popleft()
        Mkey = key_of(M)
        d = dist[Mkey]

        if d >= max_depth:
            # 打切り深さに達したら展開しない
            continue

        for i, G in enumerate(gens):
            # 即時に逆元へ戻る手（バックトラック）を禁止 → 還元語のみ列挙
            if last != -1 and i == inv_idx[last]:
                continue

            M2 = M * G
            k2 = key_of(M2)

            # 最初に単位 I に戻ったときの長さ（>0）が girth
            if k2 == Ikey and d + 1 > 0:
                if not return_word:
                    return d + 1
                # パス復元: M までの語 + 最後の i を連結
                path: List[int] = []
                cur_key = Mkey
                while True:
                    prev_key, gi = parent[cur_key]
                    if prev_key is None:
                        break
                    path.append(gi)  # 逆順で貯める
                    cur_key = prev_key
                path.reverse()
                path.append(i)
                return d + 1, path

            # 未訪問なら登録
            if k2 not in dist:
                dist[k2] = d + 1
                parent[k2] = (Mkey, i)
                q.append((M2, i))

    # 指定深さまでに閉路が見つからなかった
    return None

# ------------------------------
# 実行例
# ------------------------------
if __name__ == "__main__":
    p = 101        # 法（素数）
    n = 2          # A(n), B(n) のパラメータ
    max_depth = 32 # 打切り深さ（指数的に重くなるので注意）

    ans = cayley_girth_bfs(p, n, max_depth=max_depth, return_word=True)
    if ans is None:
        print(f"girth 未発見: p={p}, n={n}, max_depth={max_depth} まで探索")
    else:
        girth, word = ans
        print(f"girth = {girth}")
        print(f"最短閉路語（記号）: {pretty_word(word)}")

        # 検算（語の積が単位行列になることを確認）
        gens, _ = generators_A_B(n, p)
        M = apply_word(word, gens)
        print(f"検算: 語の積は単位行列か？ -> {M == identity_matrix(GF(p), 2)}")
