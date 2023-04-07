import numpy as np

def finite(q,k,t=0):
    def connect_mat(nnod,nel):
        xn = np.array(range(nnod*nnod))
        A = xn.reshape((nnod,nnod))  #nnod行,nnod列
        R = np.zeros(shape=((nnod-1)*(nnod-1),nel))   # (nnod-1)*(nnod-1)为行数,即总共的元素数，nel为列数，即4列
        for i in range((nnod-1)*(nnod-1)):
            xg = i % (nnod-1)
            yg = int(i/(nnod-1))      #生成每个元素对应的x坐标和y坐标
            a = A[xg:xg+2,yg:yg+2]  # 2*2数组，即每个元素对应的节点号
            R[i,:]=a.flatten() #修改R的第i行数字为零维的a数组
        return R
    def elemstiff2d(e,nel,h,s,connect):
        idp11 = 2/3
        idp12 = (1/6)-(1/3)
        idp13 = (1/6)-(1/3)
        idp14 = -2/6
        ke = np.diag(np.repeat(idp12,3),1)+np.diag(np.repeat(idp13,2),2)  #array是一个1维数组时，结果形成一个以一维数组为对角线元素的矩阵
        ke[0,nel-1] = idp14
        ke[1,2] = idp14
        ke = ke + ke.T + np.diag(np.repeat(idp11,nel))
        return ke
    def elemforce2d(e,nel,h,s,connect):
        fe = np.zeros(nel)
        nodes = connect[e,:]
        f = np.empty(shape=[nel,1])
        for i in range(nel):
            f[i] = df[int(nodes[i])]
        def f1(x,y):
            return f[0][0]*(x-h)*(y-h)/(h*h)
        def f2(x,y):
            return f[1][0]*x*(y-h)/(-h*h)
        def f3(x,y):
            return f[2][0]*(x-h)*y/(-h*h)
        def f4(x,y):
            return f[3][0]*x*y/(h*h)
        def gsint(f,h):
            ans = h*h*(f(0,0)+f(0,h)+f(h,0)+f(h,h))/4
            return ans
        fe = np.empty(shape=[nel,1])
        #feT = np.empty(shape=[nel,1])
        #err = np.empty(shape=[nel,1])
        #feT[0],err[0] = integrate.dblquad(f1,0,h,0,h)
        #feT[1],err[1] = integrate.dblquad(f2,0,h,0,h)
        #feT[2],err[2] = integrate.dblquad(f3,0,h,0,h)
        #feT[3],err[3] = integrate.dblquad(f4,0,h,0,h)
        fe[0] = gsint(f1,h)
        fe[1] = gsint(f2,h)
        fe[2] = gsint(f3,h)
        fe[3] = gsint(f4,h)
        return fe   
    q=np.array(q)
    k=np.array(k)
    l = len(q)
    Q = np.zeros((l*2,l*2))
    K = np.ones((l*2,l*2))
    r = int(l/2)
    Q[r:r+l, r:r+l] = q
    K[r:r+l, r:r+l] = k

    s=Q.shape[0]-1
    nnod = s+1
    el=s*s
    L = s
    h = L/s
    nod=nnod*nnod
    nel =4     #每个元素由四个节点构成
    u_b = np.zeros(2*(s+s))
    connect = connect_mat(nnod,nel)
    ebcdof = np.unique(list(range(nnod))+list(range(nnod*s,nnod*s+nnod))+list(range(nnod,nnod*(s-1)+1,nnod))+list(range(2*nnod-1,s*nnod,nnod)))
    ebcval = u_b
    bigk = np.zeros(shape=[nod,nod])
    fext = np.zeros(shape=[nod,1])
    
    df = (Q/K).flatten()
    for e in range(el):
        ke = elemstiff2d(e,nel,h,s,connect)
        fe = elemforce2d(e,nel,h,s,connect)
        sctr = connect[e]
        for i in range(nel):
            fext[int(sctr[i])] = fext[int(sctr[i])] + fe[i]
            for j in range(nel):
                bigk[int(sctr[i]),int(sctr[j])] = bigk[int(sctr[i]),int(sctr[j])] + ke[i,j]
    for i in range(len(ebcdof)):
        n = ebcdof[i]
        bigk[n,:] = 0.0
        bigk[:,n] = 0.0
        bigk[n,n] = 1.0
        fext[n] = ebcval[i]
    u_coeff = np.dot(np.linalg.inv(bigk),fext)
    u_cal_re = u_coeff.reshape(nnod,nnod)+t
    return u_cal_re[r:r+l, r:r+l]

if __name__ == '__main__':
    q=[[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,1.,1.,1.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
,[0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]]
    k=[[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,1.,1.,1.,1.,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]
,[0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1
,0.1,0.1]]
    u = finite(q,k,20)
    print(u)
