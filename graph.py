"""
SL2(F_p) 上で A(k)=[[1,k],[0,1]] と B(k)=[[1,0],[k,1]] を生成元とする
Cayley グラフ（単位元が属する連結成分）を構築・描画するスクリプト。

主な機能:
- 行列の乗算（mod p）
- 生成元 A(k), B(k)（必要なら逆元も）
- BFS による到達可能集合（語長も記録）
- NetworkX でグラフ構築し、Matplotlib で描画
"""

import argparse
from collections import deque
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

try:
    import networkx as nx
    import matplotlib.pyplot as plt
except Exception as e:
    raise SystemExit(
        "This script requires networkx and matplotlib. Please install them (e.g., pip install networkx matplotlib)."
    )


# Matrix type: (a, b, c, d) representing [[a, b], [c, d]] in F_p
Mat = Tuple[int, int, int, int]


def mat_mul(x: Mat, y: Mat, p: int) -> Mat:
    """行列 x, y（いずれも SL2(F_p) の要素として表現）を法 p で乗算する。

    引数:
    - x, y: タプル (a, b, c, d) で [[a, b], [c, d]] を表す行列
    - p: 法（素数を想定）

    戻り値:
    - x と y の積（法 p での計算結果）を同様のタプルで返す
    """
    ax, bx, cx, dx = x
    ay, by, cy, dy = y
    return (
        (ax * ay + bx * cy) % p,
        (ax * by + bx * dy) % p,
        (cx * ay + dx * cy) % p,
        (cx * by + dx * dy) % p,
    )


def mat_str(m: Mat) -> str:
    """行列表現（タプル）を可読な文字列 "[[a, b], [c, d]]" に整形する。"""
    a, b, c, d = m
    return f"[[{a}, {b}], [{c}, {d}]]"


def generators_Ak_Bk(k: int, p: int, symmetric: bool = True) -> List[Mat]:
    """SL2(F_p) における生成元 {A(k), B(k)} を返す（必要なら逆元も含める）。

    - A(k) = [[1, k], [0, 1]], B(k) = [[1, 0], [k, 1]] を F_p 上で定義。
    - symmetric=True の場合、A(k)^{-1}, B(k)^{-1} も含め、無向 Cayley グラフに対応。

    引数:
    - k: 生成元のパラメータ（p で剰余を取る）
    - p: 法（素数を想定）
    - symmetric: True なら逆元も含める

    戻り値:
    - 行列タプルのリスト（A, B, 必要なら A^{-1}, B^{-1}）
    """
    k = k % p
    A = (1 % p, k, 0, 1 % p)
    B = (1 % p, 0, k, 1 % p)
    gens = [A, B]
    if symmetric:
        A_inv = (1 % p, (-k) % p, 0, 1 % p)
        B_inv = (1 % p, 0, (-k) % p, 1 % p)
        gens.extend([A_inv, B_inv])
    return gens


def bfs_component(gens: Sequence[Mat], p: int, max_nodes: Optional[int] = None) -> Tuple[List[Mat], Dict[Mat, int], Dict[Mat, int]]:
    """生成元 gens による半群作用で、単位元が属する連結成分を BFS で列挙する。

    単位元 I から始め、左側乗（s * g）で到達できる行列を幅優先探索で集める。

    引数:
    - gens: 生成元の行列（A, B など）。無向グラフなら逆元も含める。
    - p: 法（素数を想定）
    - max_nodes: ノード数の上限（安全のための打ち切り）。None なら制限なし。

    戻り値:
    - nodes: BFS 順の行列リスト
    - index: 行列 -> 連番インデックス の辞書
    - dist: 行列 -> 単位元からの距離（語長）
    """
    I = (1, 0, 0, 1)
    q: deque[Mat] = deque([I])
    visited: Dict[Mat, int] = {I: 0}
    order: List[Mat] = [I]
    while q:
        g = q.popleft()
        d = visited[g]
        for s in gens:
            ng = mat_mul(s, g, p)  # left multiplication: s * g
            if ng not in visited:
                visited[ng] = d + 1
                order.append(ng)
                q.append(ng)
                if max_nodes is not None and len(order) >= max_nodes:
                    q.clear()
                    break
    index = {m: i for i, m in enumerate(order)}
    return order, index, visited


def build_cayley_graph(
    p: int,
    k: int,
    directed: bool = False,
    symmetric: bool = True,
    max_nodes: Optional[int] = None,
) -> Tuple["nx.Graph", List[Mat], Dict[Mat, int], Dict[Mat, int]]:
    """生成元 {A(k), B(k)} による Cayley グラフ（連結成分）を構築する。

    既定では単位元を含む連結成分を BFS で構成する。
    - A(k), B(k) が SL2(F_p) を生成するなら、これが全体の Cayley グラフに一致。
    - そうでない場合は、単位元が属する部分群の Cayley グラフとなる。

    引数:
    - p: 法（素数）
    - k: 生成元 A(k), B(k) のパラメータ
    - directed: True なら有向グラフ（A, B のみの辺）。
    - symmetric: True なら逆元も含めて無向化（既定）。
    - max_nodes: BFS の打ち切り上限（安全策）。

    戻り値:
    - G: NetworkX の Graph/DiGraph
    - nodes: BFS 順の行列リスト
    - index: 行列 -> ノード番号
    - dist: 行列 -> 単位元からの距離
    """
    gens = generators_Ak_Bk(k, p, symmetric=symmetric)
    nodes, index, dist = bfs_component(gens, p, max_nodes=max_nodes)

    G = nx.DiGraph() if directed else nx.Graph()
    # Add nodes with attributes
    for m in nodes:
        G.add_node(index[m], mat=m, label=mat_str(m), dist=dist[m])

    # Add edges according to generators
    for g in nodes:
        gi = index[g]
        for s in gens if (directed or symmetric) else gens:
            h = mat_mul(s, g, p)
            if h in index:
                hi = index[h]
                if directed:
                    G.add_edge(gi, hi)
                else:
                    if gi != hi:
                        G.add_edge(gi, hi)
    return G, nodes, index, dist


def draw_graph(
    G: "nx.Graph",
    nodes: List[Mat],
    dist: Dict[Mat, int],
    index: Dict[Mat, int],
    layout: str = "spring",
    seed: Optional[int] = 42,
    node_size: int = 80,
    edge_alpha: float = 0.35,
    label_threshold: int = 40,
    figsize: Tuple[float, float] = (8.0, 8.0),
    title: Optional[str] = None,
    cmap = plt.cm.viridis,
):
    """NetworkX グラフを描画する補助関数。

    引数:
    - G: `build_cayley_graph` で構築したグラフ
    - nodes: BFS 順の行列リスト
    - dist: 行列 -> 単位元からの距離
    - index: 行列 -> ノード番号
    - layout: レイアウト（spring/kamada_kawai/spectral/circular）
    - seed: spring レイアウトの乱数シード
    - node_size: ノードの大きさ
    - edge_alpha: 辺の透明度
    - label_threshold: ノード数がこの値以下のときだけ行列ラベルを描画
    - figsize: 図のサイズ (幅, 高さ)
    - title: タイトル文字列
    - cmap: ノード彩色用カラーマップ（距離で彩色）
    """
    # Layout
    if layout == "spring":
        pos = nx.spring_layout(G, seed=seed)
    elif layout == "kamada_kawai":
        pos = nx.kamada_kawai_layout(G)
    elif layout == "spectral":
        pos = nx.spectral_layout(G)
    elif layout == "circular":
        pos = nx.circular_layout(G)
    else:
        pos = nx.spring_layout(G, seed=seed)

    # Colors by distance from identity
    # dist is keyed by matrix; we need per-node list in node index order
    dvals = []
    for i in range(len(nodes)):
        m = nodes[i]
        dvals.append(dist[m])

    plt.figure(figsize=figsize)
    nx.draw_networkx_nodes(G, pos, node_color=dvals, cmap=cmap, node_size=node_size)
    nx.draw_networkx_edges(G, pos, alpha=edge_alpha)
    if len(nodes) <= label_threshold:
        labels = {i: G.nodes[i]["label"] for i in G.nodes}
        nx.draw_networkx_labels(G, pos, labels=labels, font_size=8)
    if title:
        plt.title(title)
    plt.axis("off")


def main():
    """コマンドライン引数を受け取り、Cayley グラフを構築・描画するエントリポイント。"""
    parser = argparse.ArgumentParser(
        description=(
            "Draw the Cayley graph of SL2(F_p) generated by A(k)=[[1,k],[0,1]] and B(k)=[[1,0],[k,1]].\n"
            "By default, the connected component of the identity is constructed via BFS."
        )
    )
    parser.add_argument("-p", "--prime", type=int, required=True, help="Prime p (modulus)")
    parser.add_argument("-k", "--k", type=int, required=True, help="Parameter k for A(k), B(k)")
    parser.add_argument(
        "--directed",
        action="store_true",
        help="Use a directed Cayley graph (edges for A, B only). Default is undirected with symmetric generators.",
    )
    parser.add_argument(
        "--no-symmetric",
        action="store_true",
        help="Do not include inverses of generators (only A, B).",
    )
    parser.add_argument(
        "--max-nodes",
        type=int,
        default=None,
        help="Optional cap for BFS exploration size (safety).",
    )
    parser.add_argument(
        "--layout",
        type=str,
        default="spring",
        choices=["spring", "kamada_kawai", "spectral", "circular"],
        help="Graph layout algorithm",
    )
    parser.add_argument("--seed", type=int, default=42, help="Seed for spring layout")
    parser.add_argument("--node-size", type=int, default=80, help="Node size for drawing")
    parser.add_argument("--edge-alpha", type=float, default=0.35, help="Edge alpha for drawing")
    parser.add_argument(
        "--label-threshold",
        type=int,
        default=40,
        help="Max node count to show matrix labels",
    )
    parser.add_argument(
        "--figsize",
        type=float,
        nargs=2,
        default=(8.0, 8.0),
        help="Figure size width height",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="If set, save the figure to this path instead of showing",
    )

    args = parser.parse_args()

    p = args.prime
    k = args.k
    symmetric = not args.no_symmetric
    directed = args.directed

    if p < 2:
        raise SystemExit("p must be a prime >= 2")
    # Note: we do not primality-check p here; user responsibility.

    G, nodes, index, dist = build_cayley_graph(
        p=p,
        k=k,
        directed=directed,
        symmetric=symmetric,
        max_nodes=args.max_nodes,
    )

    title = f"Cayley graph of <A(k),B(k)> in SL2(F_{p}) with k={k} ({'directed' if directed else 'undirected'})"
    draw_graph(
        G,
        nodes,
        dist,
        index,
        layout=args.layout,
        seed=args.seed,
        node_size=args.node_size,
        edge_alpha=args.edge_alpha,
        label_threshold=args.label_threshold,
        figsize=(args.figsize[0], args.figsize[1]),
        title=title,
    )

    if args.output:
        plt.savefig(args.output, bbox_inches="tight", dpi=200)
    else:
        plt.show()


if __name__ == "__main__":
    main()
