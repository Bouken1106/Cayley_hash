import sympy
from sympy import symbols, Poly, GF, Matrix, gcd
from sympy.polys.domains import GF
from sympy.polys.galoistools import gf_irreducible
import random as pyrandom
from sage.all import *
import time
import matplotlib.pyplot as plt
import networkx as nx
# --- 設定 ---
p = 7
n = 1
F = GF(p)

# --- 群 G = SL(2, F) の定義 ---
G = SL(2, F)

print(G)

# --- 生成元 X0, X1 の定義 ---
X0 = matrix(F, [[1, n], 
                [0, 1]])
X1 = matrix(F, [[1, 0], 
                [n, 1]])
gens = [X0, X1]

# --- 頂点の表現関数 ---
# Group élément と Matrix の両方に対応

def rep(g):
    # GAP group element なら matrix() で Sage の Matrix に変換
    m = g.matrix() if hasattr(g, 'matrix') else g
    # m.list() で行列要素をフラットなリストとして取得
    # Python の整数型に変換してタプル化
    return tuple(int(e) for e in m.list())

# --- Cayley グラフの構築 ---
C = DiGraph()               # 無向グラフ
# C = DiGraph()           # 有向グラフにする場合

# 頂点を追加
for g in G:
    C.add_vertex(rep(g))

# 辺を追加
for g in G:
    for s in gens:
        C.add_edge(rep(g), rep(g * s))

# --- グラフの描画 ---
# DiGraph C を NetworkX 用に変換
nx_g = nx.DiGraph()
for v in C.vertices():
    nx_g.add_node(v)
for e in C.edges():
    nx_g.add_edge(e[0], e[1])

# spring_layout を使用して頂点の位置を計算
pos = nx.spring_layout(nx_g)
plt.figure(figsize=(10,10))
# node_size は頂点のサイズ（面積）をピクセル単位で指定
nx.draw(nx_g, pos, node_size=20, edge_color='black', arrowsize=10)
plt.show()
