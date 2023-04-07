import numpy as np
import sympy as sp
import math
from sympy import *
import pandas as pd
def RectangleElementStiffness(E,miu,h,node_ele,p):
    s,t=sp.symbols('s t')
    x1 = node_ele[0,0]    #注意matlab与python数据索引格式的区别
    y1 = node_ele[0,1]
    x2 = node_ele[1,0]
    y2 = node_ele[1,1]
    x3 = node_ele[2,0]
    y3 = node_ele[2,1]
    x4 = node_ele[3,0]
    y4 = node_ele[3,1]
    N1 = ((1 - s) * (1 - t)) / 4
    N2 = ((1 + s) * (1 - t)) / 4
    N3 = ((1 + s) * (1 + t)) / 4
    N4 = ((1 - s) * (1 + t)) / 4
    xc = [[x1 ,x2, x3 ,x4]]      #二维数组
    xc = np.array(xc)
    yc = [[y1 ,y2, y3, y4]]
    yc = np.array(yc)
    J_m = [[0 ,1 - t ,t - s ,s - 1],
           [t - 1, 0 ,s + 1 ,- s - t],
           [s - t, - s - 1, 0 ,t + 1],
           [1 - s ,s + t, - t - 1, 0]]
    J_m= np.array(J_m)
    J=xc@J_m@yc.T/8     #矩阵相乘
    a = (y1 * (s - 1) + y2 * (-1 - s) + y3 * (1 + s) + y4 * (1 - s)) / 4
    b = (y1 * (t - 1) + y2 * (1 - t) + y3 * (1 + t) + y4 * (-1 - t)) / 4
    c = (x1 * (t - 1) + x2 * (1 - t) + x3 * (1 + t) + x4 * (-1 - t)) / 4
    d = (x1 * (s - 1) + x2 * (-1 - s) + x3 * (1 + s) + x4 * (1 - s)) / 4
    N1s=sp.diff(N1,s);     #求导
    N1t=sp.diff(N1,t);
    N2s=sp.diff(N2,s);
    N2t=sp.diff(N2,t);
    N3s=sp.diff(N3,s);
    N3t=sp.diff(N3,t);
    N4s=sp.diff(N4,s);
    N4t=sp.diff(N4,t);
    B1 = [[a * N1s - b * N1t, 0],
          [0, c * N1t - d * N1s],
          [c * N1t - d * N1s, a * N1s - b * N1t]]
    B1 = np.array(B1)
    B2 = [[a * N2s - b * N2t, 0],
          [0, c * N2t - d * N2s],
          [c * N2t - d * N2s, a * N2s - b * N2t]]
    B2 = np.array(B2)
    B3 = [[a * N3s - b * N3t, 0],
          [0, c * N3t - d * N3s],
          [c * N3t - d * N3s, a * N3s - b * N3t]]
    B3 = np.array(B3)
    B4 = [[a * N4s - b * N4t, 0],
          [0, c * N4t - d * N4s],
          [c * N4t - d * N4s, a * N4s - b * N4t]]
    B4 = np.array(B4)
    B=np.concatenate((B1,B2,B3,B4),axis=1)       #对数组进行拼接
    B = B / J;
    B = np.array(B)
    if p == 1:
        D = [[1, miu, 0],
             [miu, 1, 0],
             [0, 0, (1 - miu) / 2]]
        D = np.array(D)
        D = D.flatten()  # 降维
        miu1 = math.pow(miu, 2)     #对浮点数求平方
        miu1=E / (1 - miu1)
        D = [i *miu1  for i in D]  # 缩放
        D = np.array(D)
        D = D.reshape((3, 3))  # 还原

    # elif p == 2:
    #     D = E / (1 + miu) / (1 - 2 * miu) * [[1 - miu, miu, 0],
    #                                          [miu, 1 - miu, 0],
    #                                          [0, 0, (1 - 2 * miu) / 2]]
    D = np.array(D)
    BD = B.T@D@B*J;
    # print(BD)
    r = np.zeros((8,8))
    r= np.array(r)
    for ii in range(0, 8):
        for jj in range(0, 8):
          r[ii,jj]= integrate(integrate(BD[ii,jj], (t, -1, 1)), (s, -1, 1));    #求二重积分
    k_ele = (h * r)     #浮点数类型
    return k_ele

def assemRectangle(k_t, k_ele, node):
    d = np.zeros(8)
    d[0: 2] = [2 * node[0] - 1, 2 * node[0]]
    d[2: 4] = [2 * node[1] - 1, 2 * node[1]]
    d[4: 6] = [2 * node[2] - 1, 2 * node[2]]
    d[6: 8] = [2 * node[3] - 1, 2 * node[3]]
    d = np.array(d)
    for ii in range(0, 8):
        for jj in range(0, 8):
            k_t[int(d[ii] - 1), int(d[jj] - 1)] = k_t[int(d[ii] - 1), int(d[jj] - 1)] + k_ele[ii, jj]
    return k_t

def RectangleElementStress(E,miu,node_ele,u1,p):
    s, t = sp.symbols('s t')
    x1 = node_ele[0, 0]  # 注意matlab与python数据索引格式的区别
    y1 = node_ele[0, 1]
    x2 = node_ele[1, 0]
    y2 = node_ele[1, 1]
    x3 = node_ele[2, 0]
    y3 = node_ele[2, 1]
    x4 = node_ele[3, 0]
    y4 = node_ele[3, 1]
    N1 = ((1 - s) * (1 - t)) / 4
    N2 = ((1 + s) * (1 - t)) / 4
    N3 = ((1 + s) * (1 + t)) / 4
    N4 = ((1 - s) * (1 + t)) / 4
    xc = [x1, x2, x3, x4]
    yc = [y1, y2, y3, y4]
    J_m = [[0, 1 - t, t - s, s - 1],
           [t - 1, 0, s + 1, - s - t],
           [s - t, - s - 1, 0, t + 1],
           [1 - s, s + t, - t - 1, 0]]
    J_m= np.array(J_m)
    xc = np.array(xc)
    yc = np.array(yc)
    J = xc @ J_m @ yc.T / 8

    a = (y1 * (s - 1) + y2 * (-1 - s) + y3 * (1 + s) + y4 * (1 - s)) / 4
    b = (y1 * (t - 1) + y2 * (1 - t) + y3 * (1 + t) + y4 * (-1 - t)) / 4
    c = (x1 * (t - 1) + x2 * (1 - t) + x3 * (1 + t) + x4 * (-1 - t)) / 4
    d = (x1 * (s - 1) + x2 * (-1 - s) + x3 * (1 + s) + x4 * (1 - s)) / 4
    N1s = sp.diff(N1, s);  # 求导
    N1t = sp.diff(N1, t);
    N2s = sp.diff(N2, s);
    N2t = sp.diff(N2, t);
    N3s = sp.diff(N3, s);
    N3t = sp.diff(N3, t);
    N4s = sp.diff(N4, s);
    N4t = sp.diff(N4, t);
    B1 = [[a * N1s - b * N1t, 0],
          [0, c * N1t - d * N1s],
          [c * N1t - d * N1s, a * N1s - b * N1t]]
    B1 = np.array(B1)
    B2 = [[a * N2s - b * N2t, 0],
          [0, c * N2t - d * N2s],
          [c * N2t - d * N2s, a * N2s - b * N2t]]
    B2 = np.array(B2)
    B3 = [[a * N3s - b * N3t, 0],
          [0, c * N3t - d * N3s],
          [c * N3t - d * N3s, a * N3s - b * N3t]]
    B3 = np.array(B3)
    B4 = [[a * N4s - b * N4t, 0],
          [0, c * N4t - d * N4s],
          [c * N4t - d * N4s, a * N4s - b * N4t]]
    B4 = np.array(B4)
    B = np.concatenate((B1, B2, B3, B4), axis=1)        #拼接矩阵
    B = B / J;
    B = np.array(B)
    if p == 1:
        D = [[1, miu, 0],
             [miu, 1, 0],
             [0, 0, (1 - miu) / 2]]
        D=np.array(D)
        D = D.flatten()     #降维
        miu1 = math.pow(miu, 2)
        miu1 = E / (1 - miu1)
        D = [i * miu1 for i in D]      #缩放
        print(D)
        D = np.array(D)
        D = D.reshape((3, 3))       #还原
    # elif p == 2:
    #     D = E / (1 + miu) / (1 - 2 * miu) * [[1 - miu, miu, 0],
    #                                          [miu, 1 - miu, 0],
    #                                          [0, 0, (1 - 2 * miu) / 2]]
    D = np.array(D)
    str = D @ B @ u1;  # 单元内的应力场
    print(str)
    wcent=np.zeros((3,1))
    # print(str[0,0])
    # str1= (-26.6666666666669*s - 800.000000000002*t - 373.337156401729).subs(s,0).subs(t, 0);
    # print(str1)
    for i in range(3):
        wcent[i, 0] = str[i, 0].subs(s, 0);  # 单元中心应力值
        wcent[i, 0] = str[i, 0].subs(t, 0);
        # w = wcent;   #浮点数类型
    return wcent


# # a=np.arange(9)
# # print(a)
# file = np.load('C:/Users/sexylover/Downloads/a.npy')
#
# np.set_printoptions(threshold=np.inf)    #完全显示矩阵的所有值   [z][y][x]
# # print(file)
# # c=file.shape              #显示矩阵的形状
# # print(c)
# # np.savetxt('C:/Users/sexylover/Downloads/a.txt')     #将ppy文件保存为txt文件
# example1=file[:,12,:]     #单元受力情况
# hangshu=len(example1)
# nnod=hangshu+1
# # print(file[:,12,:])
# # print(file[:,12,:].shape)
# # nnod = 6
# # nel = 4
# # xn = np.array(range(nnod * nnod))
# # A = xn.reshape((nnod, nnod))  # nnod行,nnod列
# # R = np.zeros(shape=((nnod - 1) * (nnod - 1), nel))  # (nnod-1)*(nnod-1)为行数,即总共的单元数，nel为列数，即4列
# # for i in range((nnod - 1) * (nnod - 1)):
# #     xg = i % (nnod - 1)
# #     yg = int(i / (nnod - 1))  # 生成每个单元对应的x坐标和y坐标
# #     a = A[xg:xg + 2, yg:yg + 2]  # 2*2数组，即每个单元对应的节点号
# #
# #     R[i, :] = a.flatten()      #降维
# # print(R)
# A = np.zeros(shape=(nnod, nnod ))
# shoulidianshu=0
# for i in range(hangshu * hangshu):
#     xg = i % hangshu
#     yg = int(i / hangshu)  # 生成每个单元对应的x坐标和y坐标
#     # a = A[xg:xg + 2, yg:yg + 2]  # 2*2数组，即每个单元对应的节点号
#     if example1[xg,yg]!=0:
#         A[xg+1,yg+1]=example1[xg,yg]
#         shoulidianshu=shoulidianshu+1
# print(A)     #节点受力情况
#
# R = np.zeros(shape=(nnod * nnod , 4))
# for i in range(nnod * nnod):
#     xg = i % nnod
#     yg = int(i / nnod)
#     R[i,0]=i
#     R[i, 1] = xg
#     R[i, 2] = yg
#     R[i, 3] = 0
# print(R)      #节点信息，第一列为节点编号，2~4列分别为x,y,z方向坐标
#
# xn = np.array(range(nnod * nnod))
# B = xn.reshape((nnod, nnod))
# Z= np.zeros(shape=(shoulidianshu, 5))
# count=0
# for i in range(nnod * nnod):
#     xg = i % nnod
#     yg = int(i / nnod)  # 生成每个节点对应的x坐标和y坐标
#     if A[xg,yg]!=0:
#         Z[count,0]=count   #参与计算的单元编号
#         row = np.array([xg-1,xg-1,xg,xg])
#         col =np.array([yg-1,yg,yg,yg-1])
#         Z[count, 1:5] = B[row,col]       #取特定列，特定行的节点
#         count = count + 1
# print(Z)     #单元信息，第一列为单元编号，后面各列为单元上的节点号码
node = [[1, 0, 0, 0],
        [2, 1, 0, 0],
        [3, 2, 0, 0],
        [4, 0, 1, 0],
        [5, 1, 1, 0],
        [6, 2, 1, 0]]
node= np.array(node)
ele=[[1,1,2,5,4],
     [2,2,3,6,5]]
ele= np.array(ele)
E=2.1e4;             #弹性模量
t=0.025;               #单元厚度
miu=0.3;          #泊松比
num_ele = np.size(ele, 0); # 单元数
num_node = np.size(node, 0); # 节点数
dof=num_node *2;       #自由度数，梁单元每个节点有2个自由度，x,y方向位移
f=np.ones((dof,1))*1e8;             #结构整体外载荷矩阵，整体坐标系下
f_loc=np.zeros((8,1));              #单元外载荷矩阵，局部坐标系下
u=np.ones((dof,1))*1e6;             #位移矩阵
K=np.zeros((dof,dof));                  #总体刚度矩阵
for i in range(num_ele):
    k_ele=RectangleElementStiffness(E,miu,t,node[ele[i,1:5]-1,1:3],1);
    K=assemRectangle(K,k_ele,ele[i,1:5]);
f[10]=3
f[11]=4
for i in range(2*num_node):
    if f[i] == 1e8:
        f[i] = 0
f[0]=1e8
f[1]=1e8
f[6]=1e8
f[7]=1e8

u[0]=0;
u[1]=0;
u[6]=0;
u[7]=0;

index=[];      #未知自由度的索引
p=[];          #未知自由度对应的节点力矩阵

for i in range(dof):
    if u[i]!=0:
        index=index+[i]
        p = np.concatenate((p, f[i]), axis=0)
p=np.array([p])
p=p.T

print(index)
index=np.array(index)
print(p)
# for i in range(dof):
K1=K[index, :]
K1=K1[:, index]
print(K1)
u[index] =np.dot(np.linalg.inv(K1), p)
f=K@u;
stress=np.zeros((num_ele,3));
print(2*ele[1,1])
for i in range(num_ele):
    u1=np.zeros((8,1))
    u1[0]=u[2*ele[i,1]-2]
    u1[1]=u[2*ele[i,1]]
    u1[2]=u[2*ele[i,2]-2]
    u1[3]=u[2*ele[i,2]]
    u1[4]=u[2*ele[i,3]-2]
    u1[5]=u[2*ele[i,3]]
    u1[6]=u[2*ele[i,4]-2]
    u1[7]=u[2*ele[i,4]]
    u1= np.array(u1)
    stress[i,:]=np.array(RectangleElementStress(E,miu,node[ele[i,1:5]-1,1:3],u1,1)).T;
print(stress)