from math import cos, sin
from os import P_NOWAIT
from tkinter import N
import numpy as np

class Node:
    # 节点类
    def __init__(self, type, P = 0, Q = 0, U = 1, delta = 0) :
        self.type = type
        self.P = P
        self.Q = Q
        self.U = U
        self.delta = delta
        self.print_node()
    
    def print_node(self):
        print("*********************************************")
        print("当前节点类型为: " + self.type + " 型")
        print("当前节点P= " + str(self.P))
        print("当前节点Q= " + str(self.Q))
        print("当前节点U= " + str(self.U))
        print("当前节点δ= " + str(self.delta))
        print("*********************************************")

def DeltaPQ():
    # 功率误差
    deltaP = np.zeros([n,1])
    deltaQ = np.zeros([n,1])
    for i in range(n):
        sumP = 0
        sumQ = 0
        for j in range(n):
            sumP += Un[j]*(
                G[i,j] * cos(deltaN[i,j]) + 
                B[i,j] * sin(deltaN[i,j]))
            sumQ += Un[j]*(
                G[i,j] * sin(deltaN[i,j]) - 
                B[i,j] * cos(deltaN[i,j]))
        deltaP[i,0] = Pn[i] - Un[i]*sumP
        deltaQ[i,0] = Qn[i] - Un[i]*sumQ

    return np.append(deltaP,deltaQ,axis=0)   

def Jacoby():
    # 初始化分块雅可比矩阵
    H = np.zeros([n,n])
    M = np.zeros([n,n])
    N = np.zeros([n,n])
    L = np.zeros([n,n])
    
    for i in range(n):
        for j in range(n):
            if i != j :
                # 非对角
                H[i,j] = Un[i] * Un[j] * (
                    G[i,j] * sin(deltaN[i,j]) - 
                    B[i,j] * cos(deltaN[i,j]))
                M[i,j] = - Un[i] * Un[j] * (
                    G[i,j] * cos(deltaN[i,j]) + 
                    B[i,j] * sin(deltaN[i,j]))
                N[i,j] = -M[i,j]
                L[i,j] = H[i,j]
            else :
                # 对角
                H[i,j] = - Un[i]*Un[i]*B[i,j] - Qn[i]
                N[i,j] = + Un[i]*Un[i]*G[i,j] + Pn[i]
                M[i,j] = - Un[i]*Un[i]*G[i,j] + Pn[i]
                L[i,j] = - Un[i]*Un[i]*B[i,j] + Qn[i]
        
                
    #返回拼接雅可比矩阵
    return np.append(
        np.append(H,N,axis=0),
        np.append(M,L,axis=0),
        axis=1)     

def DeltaUdot(pq,j):
    # 电压幅值、相位误差
    return np.matmul(np.linalg.inv(j),pq)

print("本例采用牛-拉法机算3节点潮流一次迭代")
print("=====================================================")

nodes = [
    Node("BN",U = 1, delta = 0),
    Node("PV",P = 0.7, U = 1.05),
    Node("PQ",P = -2.8, Q = -1.2)]
n = 3
deltaN = np.zeros([n,n])
Pn = np.array([node.P for node in nodes])
Qn = np.array([node.Q for node in nodes])

# 各节点间支路阻抗
y = np.array([
    0, 1/0.1j, 1/0.1j,
    1/0.1j, 0, 1/0.1j,
    1/0.1j, 1/0.1j, 0]).reshape([n,n])
# print(Z.shape)
# 由于我的矩阵法电路知识已经不足以支撑我列出ZY之间的矩阵关系式
# 所以我选择直接循环计算
Y = np.zeros([n,n]).astype("complex")
for i in range(n):
    for k in range(n):
        if i == k :
            #Yii
            Y[i,k] = np.sum(y[i,:])
        else :
            #Yik
            Y[i,k] = -y[i,k]
print("=====================================================")
print("节点导纳矩阵Y = ")
print(Y)
print("=====================================================")
# 实虚分离构造GB矩阵便于计算
# 注意B矩阵为虚部!(即不含j)
G = np.real(Y)
B = np.imag(Y)

# 迭代次数
iterations = 1
# 容许误差
epsilon = 1

for k in range(iterations):
    
    Un = np.array([node.U for node in nodes])
    for i in range(n):
        for k in range(n):
            deltaN[i,k] = nodes[i].delta - nodes[k].delta
            
    deltaPQ = DeltaPQ()
    print("====================================================")
    print("当前功率误差向量 = ")
    print(deltaPQ)
    print("====================================================")
    if max(abs(deltaPQ.copy())) < epsilon:
        # 收敛判据
        break
    J = Jacoby()
    print("====================================================")
    print("当前雅各比矩阵J = ")
    print(J)
    print("====================================================")
    deltaPQ = np.delete(deltaPQ , [0,3,4] , axis=0 )
    J = np.delete(J , [0,3,4] , axis=0 )
    J = np.delete(J , [0,3,4] , axis=1 )
    correctionU = DeltaUdot(deltaPQ.copy(), J.copy())
    ("====================================================")
    print("整理后功率误差向量 = ")
    print(deltaPQ)
    print("====================================================")
    ("====================================================")
    print("整理后雅各比矩阵J = ")
    print(J)
    print("====================================================")
    ("====================================================")
    print("电压误差矩阵 = ")
    print(correctionU)
    print("====================================================")
    for node in nodes:
        if node.type == "BN":
            # 平衡节点，跳过
            continue
        elif node.type == "PV":
            # PV节点，仅修改δ
            node.delta += correctionU[0,0]
        else :
            # PQ节点，修改U、δ
            node.U *= 1+correctionU[2,0]
            node.delta += correctionU[1,0]
            
print(np.array([node.U for node in nodes]))
print(np.array([node.delta for node in nodes]))