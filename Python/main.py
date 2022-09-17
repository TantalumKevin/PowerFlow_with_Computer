import numpy as np

class node:
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

       
            
print("本例采用牛-拉法机算3节点潮流一次迭代")
print("=====================================================")

nodes = [
    node("BN",U = 1, delta = 0),
    node("PV",P = 0.7, U = 1.05),
    node("PQ",P = -2.8, Q = -1.2)]

Z = np.array([
    0, 0.1j, 0.1j,
    0.1j, 0, 0.1j,
    0.1j, 0.1j, 0]).reshape([3,3])

#print(Z.shape)
