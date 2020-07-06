import numpy as np
import matplotlib.pyplot as plt
import cv2
import scipy.io as sio
import time as t
import networkx as nx

from multiprocessing import Pool


    level = [(0,-1)]
    marked = [(0,-1)]
    C = [0]
    while level != []:
        next_level = []
        for u in level:
            for i in range(len(Gf[0])):
                if( Gf[u[0]][i] != 0 and i not in C ):
                    if( i == len(Gf[0])-1 ):
                        def findAugmPathBreadth(Gf):
                        P = [i]
                        aux = u[0]
                        while( aux != -1 ):
                            P.append(aux)
                            for c in marked:
                                if(aux == c[0]):
                                    aux = c[1]
                                    marked.remove(c)
                        return [P[len(P)-i-1] for i in range(len(P))]
                    marked.append((i,u[0]))
                    C.append(i)
                    next_level.append((i,u[0]))
        level = [next_level[i] for i in range(len(next_level))]
    return []

def hungarianMethodStep1(A):
    for i in range(len(A[0])):
        aux = min(A[i])
        A[i] -= aux
    return A

def hungarianMethodStep2(A):
    At = np.transpose(A)
    for i in range(len(A[0])):
        aux = min(At[i])
        for j in range(len(A[0])): A[j][i] -= aux
    return A

def hungarianMethodStep3(A):
    B = A.copy()
    n = len(A[0])
    M = []
    Z =np.zeros((2*n+2,2*n+2))
    for i in range(n):
        for j in range(n):
            if(B[i][j] == 0): Z[1+i][n+j+1] = 1
    #inicializa aristas de s y t
    for i in range(n):
        Z[0][1+i] = 1
        Z[n+1+i][2*n+1] = 1
    P = findAugmPathBreadth(Z)
    while(P != []):
        for j in range(len(P)-1):
            Z[P[j]][P[j+1]] -= 1
            Z[P[j+1]][P[j]] += 1
        P = findAugmPathBreadth(Z)
    mod = 0
    for i in range(n,2*n+1):
        for j in range(1,n+1):
            if(Z[i][j] != 0):
                M.append((i-n,j))
                mod += 1
    if(mod == n): return M
    columns = []
    rows = []
    for i in range(1,n+1):
        if(i not in [M[j][1] for j in range(len(M))]): rows.append(i)
    for k in rows:
        for i in range(n):
            if(B[k-1][i] == 0 and i+1 not in columns):
                columns.append(i+1)
                for j in range(len(M)):
                    if(M[j][0] == i+1 and M[j][1] not in rows): rows.append(M[j][1])
    Left = []
    for i in range(1,n+1):
        for j in range(1,n+1):
            if(j not in columns and i in rows): Left.append((B[i-1][j-1],i-1,j-1))
    val = min([Left[i][0] for i in range(len(Left))])
    for K in Left: B[K[1]][K[2]] -= val

    M1 = hungarianMethodStep3(B)
    return M1

def hungarianAlgorithm(A):
    hungarianMethodStep1(A)
    hungarianMethodStep2(A)
    return hungarianMethodStep3(A)

def norm(p,q):
    N = ( (p[0]-q[0])**2 + (p[1]-q[1])**2 )**(0.5)
    #print(N)
    return N

def histograma(i,P):
    radBins = 8        # Cantidad de cajas radiales (anillos)
    angleBins = 15      # Cantidad de cajas angulares (quesitos)
    PA = 1.0*(np.copy(P))
    """
    X =[0,0]
    for i in range(len(PA)):
        X[0] += 1.0*PA[i][0]/len(PA)
        X[1] += 1.0*PA[i][1]/len(PA)
    """
    # Normaliza la forma con centro en X
    Z = P[i]
    #print(PA)
    Max = max([norm(P[j],Z) for j in range(len(P))])
    for j in range(len(PA)):
        a = float(P[j][0]-Z[0])
        b = float(P[j][1]-Z[1])
        PA[j][0] = a/(Max+1)
        PA[j][1] = b/(Max+1)
    # Inicializa el histograma
    H = np.zeros([radBins,angleBins])

    # Llena los bins (matriz de presencia)
    for j in range(len(PA)):
        if(j != i):
            theta = np.arctan2(PA[j][1],PA[j][0])
            r = norm(PA[i],PA[j])*(2**(radBins))
            #print(theta)
            if(r != 0):
                if np.log2(r) < 0:
                    H[0][int(1.0*angleBins*((theta+np.pi)%(2*np.pi))/(2*np.pi))] += 1.0
                else:
                    H[int(np.log2(r))][int(1.0*angleBins*((theta+np.pi)%(2*np.pi))/(2*np.pi))] += 1.0
    return H/float(sum(sum(H)))

def costFunction(ptInFirstShape,ptInSecondShape,H):
    #print(shape2)

    H1 = H[0][ptInFirstShape]
    H2 = H[1][ptInSecondShape]

    C = 0
    for k1 in range(len(H1)):
        for k2 in range(len(H1[0])):
            if H1[k1][k2]+H2[k1][k2] != 0:
                C += ((H1[k1][k2]-H2[k1][k2])**2)/(H1[k1][k2]+H2[k1][k2])
    return C/2.0

def calculaHistogramas(X,Y):
    size = len(X)
    mul8 = size//8
    res = size-mul8*8
    H = [[],[]]
    for i in range(mul8):
        with Pool(8) as p:
            h0 = p.starmap(histograma,[(8*i,X),(8*i+1,X),(8*i+2,X),(8*i+3,X),(8*i+4,X),(8*i+5,X),(8*i+6,X),(8*i+7,X)])
            h1 = p.starmap(histograma,[(8*i,Y),(8*i+1,Y),(8*i+2,Y),(8*i+3,Y),(8*i+4,Y),(8*i+5,Y),(8*i+6,Y),(8*i+7,Y)])
        for j in range(8):
            H[0].append(h0[j])
            H[1].append(h1[j])
    for i in range(mul8,size):
        H[0].append(histograma(i,X))
        H[1].append(histograma(i,Y))
    return H

def shapeContext(X, Y):
    size = len(X)
    C = np.zeros([size,size])
    #print(X, Y)

    H = calculaHistogramas(X,Y)

    for i in range(size):
        for j in range(size):
            C[i][j] = costFunction(i,j,H)
    return C

def shapenorm(X,Y):
    sizeX = len(X)
    sizeY = len(Y)
    size = 0
    if sizeX > sizeY:
        for _ in range(sizeX-sizeY):
            Y.append(Y[0])
        size = sizeX
    elif sizeY > sizeX:
        for _ in range(sizeY-sizeX):
            X.append(X[0])
        size = sizeY
    else: size = sizeX


    C = shapeContext(X,Y)
    A = hungarianAlgorithm(C)
    norm = 0
    for (i,j) in A:
        norm += C[i-1][j-1]
    return [A,norm]



scale = 9

for j in range(6,262):
    Normas = []
    cosa1 = np.array(sio.loadmat("data/"+str(j)+".mat")['x'])
    for k in range(j+1,262):
        cosa2 = np.array(sio.loadmat("data/"+str(k)+".mat")['x'])

        Dis1 = len(cosa1)
        Dis2 = len(cosa2)

        scale1 = scale
        scale2 = scale
        if Dis1 > Dis2:
            for i in range(scale):
                if Dis1//scale <= Dis2//(scale-i):
                    break
                else: scale2 = (scale-i)
        elif Dis1 < Dis2:
            for i in range(scale):
                if Dis2//scale <= Dis1//(scale-i):
                    break
                else: scale1 = (scale-i)
        print(Dis1//scale1,Dis2//scale2)

        t1 = t.time()
        Res = shapenorm([[int(cosa1[i][0]),int(cosa1[i][1])] for i in range(0,len(cosa1),scale1)],[[int(cosa2[i][0]),int(cosa2[i][1])] for i in range(0,len(cosa2),scale2)])
        t2 = t.time()


        A = Res[0]
        dist = Res[1]
        Normas.append(dist)

        print('La distancia entre las imágenes '+str(j)+' y '+str(k)+' es: ', dist)
        print("Se tardó ", t2-t1, " segundos")
    file = open('distancias'+str(j)+'.txt',"w+")
    for k in range(len(Normas)):
        file.write("D("+str(j)+","+str(k+3)+") = "+str(Normas[k])+"\n")
    file.close()

    print(Normas)
    t.sleep(10)

"""
cosa1 = np.array(sio.loadmat("data/2.mat")['x'])
cosa2 = np.array(sio.loadmat("data/7.mat")['x'])

Dis1 = len(cosa1)
Dis2 = len(cosa2)

scale1 = scale
scale2 = scale
if Dis1 > Dis2:
    for i in range(scale):
        if Dis1//scale <= Dis2//(scale-i):
            break
        else: scale2 = (scale-i)
elif Dis1 < Dis2:
    for i in range(scale):
        if Dis2//scale <= Dis1//(scale-i):
            break
        else: scale1 = (scale-i)
print(Dis1//scale1,Dis2//scale2)

t1 = t.time()
Res = shapenorm([[int(cosa1[i][0]),int(cosa1[i][1])] for i in range(0,len(cosa1),scale1)],[[int(cosa2[i][0]),int(cosa2[i][1])] for i in range(0,len(cosa2),scale2)])
t2 = t.time()


A = Res[0]
dist = Res[1]
Pl = []

print('La distancia entre las imágenes es: ', dist)
print("Se tardó", t2-t1)



for i in range(len(A)):
    try:
        Pl.append((cosa1[scale*A[i][1]-scale][0],cosa1[scale*A[i][1]-scale][1]))
    except: Pl.append((cosa1[0][0],cosa1[0][1]))
    try:
        Pl.append((cosa2[scale*A[i][0]-scale][0],cosa2[scale*A[i][0]-scale][1]))
    except: Pl.append((cosa2[0][0],cosa2[0][1]))
for i in range(0,len(Pl),2):
    plt.subplot(1,1,1),plt.plot([Pl[i][0],Pl[i+1][0]],[Pl[i][1],Pl[i+1][1]], '--')

plt.subplot(1,1,1),plt.plot(cosa1[:,0],cosa1[:,1],'-b')
plt.subplot(1,1,1),plt.plot(cosa2[:,0],cosa2[:,1],'-r')

plt.show()
"""
