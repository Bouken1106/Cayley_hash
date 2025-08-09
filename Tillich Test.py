import sympy
from sympy import symbols, Poly, Matrix, gcd
from sympy.polys.domains import GF
from sympy.polys.galoistools import gf_irreducible
import random as pyrandom
from sage.all import *






# GF(2)上の次数3既約多項式を1つランダムに生成
def Random_Irreducible_Polynomial(degree):
    F = GF(2)                          
    R = PolynomialRing(F, 'x')         
    x_poly = R.gen()                     
    while True:
        f_candidate = R.random_element(degree)
        if f_candidate.is_irreducible():
            f = f_candidate
            break
    return f


F = GF(2)                          
R = PolynomialRing(F, 'x')         
x_poly = R.gen()                     
while True:
    f_candidate = R.random_element(degree=3)
    if f_candidate.is_irreducible():
        f = f_candidate
        break
print(f)
print(Random_Irreducible_Polynomial(3))  # 例えば x^3 + x + 1 など
