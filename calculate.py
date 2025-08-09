import sympy
from sympy import symbols, Poly, GF, Matrix, gcd, simplify, cancel
from sympy.polys.domains import GF
from sympy.polys.galoistools import gf_irreducible
import random as pyrandom
from sage.all import *
import time
import matplotlib.pyplot as plt
import networkx as nx

A = Matrix(2,[[1,2],
              [0,1]])
B = Matrix(2,[[1,0],
              [2,1]])

x = symbols('x')
p = 3
n = 4
d = 2
expr = (x**(p**n) - x)/(x**(p**d) - x)

# 方法1: simplify を利用
result1 = simplify(expr)
print("simplify:", result1)

# 方法2: cancel を利用（分数の約分を行う）
result2 = cancel(expr)
print("cancel:", result2)

print(sympy.__file__)
