from sage.all import *
import time
import random
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')
import numpy as np
import sympy
from sympy import primerange, nextprime
import sys
sys.setrecursionlimit(10000)  # デフォルトは1000

def gcd(a, b):
    a, b = int(a), int(b)
    while b:
        a, b = b, a % b
    return a

def Euclid(a,b,que): #a,bは互いに素,aは0でないとする
    if b == 0:
        return "".join(que)
    if b == 1:
        que.extend(["x"] * (a - 1)) #x = matrix([[1,1],[0,1]]), y = matrix([[1,0],[1,1]])
        que.append("y")
        return "".join(que)
    if a >= b:
        r = a % b
        q = a // b
        que.extend(["x"] * q)
        return Euclid(r,b,que)
    if b >= a:
        r = b % a
        q = b // a
        que.extend(["y"] * q)
        return Euclid(a,r,que)

    
#ランダムにSL(2,Z+)の元を生成

def random_SL2Z_plus(word_length):
    X = Matrix([[1, 1], [0, 1]])
    Y = Matrix([[1, 0], [1, 1]])
    M = Matrix([[1, 0], [0, 1]])
    for _ in range(word_length):
        mat = random.choice([X, Y])
        M = M * mat
    return M

print(random_SL2Z_plus(10))

#Euclidの互除法を用いた攻撃方法


def decrypt_Euclid(M, que):  # 再帰を除去したバージョン
    while True:
        if M[1][0] == 0:
            que.extend(["X"] * M[0][1])
            return "".join(que)
        if M[1][0] == 1:
            q = M[0][0] // M[1][0]
            M[0][0] = 1
            M[0][1] = M[0][1] - (q - 1) * M[1][1]
            que.extend(["X"] * (q - 1))
            if M[0][1] == 0:
                que.append("Y")
            else:
                que.append("Y")
                que.extend(["X"] * (M[0][1]))
            return "".join(que)
        if M[0][0] >= M[1][0]:
            q = M[0][0] // M[1][0]
            M[0][0] = M[0][0] - q * M[1][0]
            M[0][1] = M[0][1] - q * M[1][1]
            que.extend(["X"] * q)
        else:  # M[1][0] >= M[0][0]
            q = M[1][0] // M[0][0]
            M[1][0] = M[1][0] - q * M[0][0]
            M[1][1] = M[1][1] - q * M[0][1]
            que.extend(["Y"] * q)

print(decrypt_Euclid([[1,1],[0,1]],[]))

#総当たりを用いた攻撃方法
def decrypt_brute_force(M):
    X = Matrix([[1,1],[0,1]])
    Y = Matrix([[1,0],[1,1]])
    if M == Matrix([[1,0],[0,1]]):
        return "Identity matrix"
    for i in range(2**1,2**256):
        E = Matrix([[1,0],[0,1]])
        binary = bin(i)[3:]
        que = []
        for j in range(len(binary)):
            if binary[j] == "1":
                E = E * X
                que.append("X")
            else: 
                E = E * Y
                que.append("Y")
        if E == M:
            return "".join(que)
    return "intractable to factor"

print(decrypt_brute_force(Matrix([[33,23],[43,30]])))

#argumentsの長さのword lengthに対して、ユークリッドの互除法を用いてどれだけ時間がかかるかを100回行い平均をとる
arguments = [0,1000,2000,3000]
time_list = []
for i in arguments:
    x = 0
    for k in range(100):
        M = random_SL2Z_plus(i)
        M = [list(M[0]),list(M[1])]
        
        start = time.perf_counter()
        decrypt_Euclid(M,[])
        end = time.perf_counter()
        x = x + end - start
    time_list.append(x/100)

x_axis = np.arange(0,3000,1000)
plt.xlabel("word size")
plt.ylabel("time(s)")
plt.plot(arguments, time_list)
plt.show()
plt.close()  # プロットをクリア

