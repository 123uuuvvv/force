import numpy as np
import math
import sympy as sp
from sympy import *
#
# # # d=np.zeros(100)
# # # e=d.reshape((2,50))
# # # a=range(0,100)
# # # d[0:100]=a;
# # # c=len(d)
# # # print(c,d)
# # # node = [[1, 0, 0, 0],
# # #         [2, 1, 0, 0],
# # #         [3, 2, 0, 0],
# # #         [4, 0, 1, 0],
# # #         [5, 1, 1, 0],
# # #         [6, 2, 1, 0]]
# # # node= np.array(node)
# # # ele=[[1,1,2,5,4],
# # #      [2,2,3,6,5]]
# # # ele= np.array(ele)
# # # c=node[ele[0,1:5]-1,1:3]
# # # # print(c)
# # # d=node.T
# # # # print(d)
# # # xc=[[1, 1 ,1, 1]];
# # # yc=[[2, 2, 2 ,2]];
# # # J_m=[[0, 1 ,0, 1],
# # #      [4,0 ,5 ,7],
# # #      [5 ,3, 0 ,5],
# # #      [3 ,2 ,1 ,0]]
# # # xc = np.array(xc)
# # # yc = np.array(yc)
# # # J_m= np.array(J_m)
# # # J=xc@J_m@yc.T/8
# # # print(J)
# # B1=[[1 ,0],
# #    [ 0 ,1],
# #     [1, 1]]
# # B1= np.array(B1)
# # B1=B1.flatten()
# # E=2.1e4;
# # t=0.025;
# # miu=0.3;
# # miu1=math.pow(0.3,2)
# # miu1=E / (1 - miu1)
# # # B2=[1,2,3]
# # # B=np.concatenate((B1,B1,B1,B1),axis=1)
# # # B = B / J;
# # # print(B)
# # # B1=[i*miu1 for i in B1]
# # # B1= np.array(B1)
# # # B1= B1.reshape((3, 2))
# #
# # print(miu1)
# s, t = sp.symbols('s t')
# BD = [[s + t, t],
#       [t, s]]
# BD = np.array(BD)
# r = np.zeros((2,2))
# r = np.array(r)
# for ii in range(0, 2):
#     for jj in range(0, 2):
#         r[ii, jj] = integrate(integrate(BD[ii, jj], (t, -2, 1)), (s, -2, 1));
# r = np.array(r)
# r=2*r
# print(r)
# node=[1,2,5,4]
# d = np.zeros(8)
# d[0: 2] = [2 * node[0] - 1, 2 * node[0]]
# d[2: 4] = [2 * node[1] - 1, 2 * node[1]]
# d[4: 6] = [2 * node[2] - 1, 2 * node[2]]
# d[6: 8] = [2 * node[3] - 1, 2 * node[3]]
# print(d)
# index=[];      #未知自由度的索引
# p=[];          #未知自由度对应的节点力矩阵
# print(type(index))
# u=[4,6,8,9,7,8,0,7]
# f=[4,6,8,9,7,8,0,7]
#
# for i in range(8):
#     if u[i]!=0:
#         p = np.concatenate((p, [f[i]]), axis=0)
#
#         # p = np.concatenate((p, i), axis=0)
# p=np.array([p])
# p=p.T
# print(p)
a=[0:5]
print(a)