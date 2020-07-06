import numpy as np
import matplotlib.pyplot as plt
import cv2

X = 256
Y = 256
D = [2,2]
N = 10000





epsilon = 0.1
scale = 1
Th = np.pi
Mat=[[np.cos(Th),-np.sin(Th)],[np.sin(Th),np.cos(Th)]]

def findAdj(A,T,v): # Find adjacent nodes of node v contained in A
    adj = []        # excluding the nodes in the subtree T
    V = [T[i][2] for i in range(len(T)) ]
    for i in range(A[0].size):
        if( A[v[2]][i] != float('inf') and i not in V ): adj.append([A[v[2]][i],(v[2],i)])
    return adj
def getPath(P, nod0, nodn, PHI, ind = -1):
    Path = []
    if( nod0 == nodn ): return [nod0]
    if( PHI == 'I' ):
        if( P[nod0][nodn] == -1 ):
            if( ind == -1 ): return [nod0,nodn]
            elif( ind == 1 ): return [nodn]
        elif( P[nod0][nodn] != nod0 and P[nod0][nodn] != nodn ):
            Path += getPath(P,nod0,P[nod0][nodn],PHI,-1)
            Path += getPath(P,P[nod0][nodn],nodn,PHI,1)
            return Path
        else:
            return [nod0, nodn]
    elif( PHI == 'P' ):
        if( P[nod0][nodn] == nod0 ):
            Path = [int(nod0),int(nodn)]
            return Path
        else:
            Path = getPath(P,nod0,P[nod0][nodn],PHI)
            Path.append(int(nodn))
            return Path
    elif( PHI == 'H' ):
        if( P[nod0][nodn] == nodn ):
            Path = [int(nod0),int(nodn)]
            return Path
        else:
            Path = getPath(P,P[nod0][nodn],nodn,PHI)
            Path.insert(0,int(nod0))
            return Path

def findAugmPathDepth(Gf, k=0, P = [0]):
    for i in range(len(Gf[0])):
        if( Gf[k][i] != float('inf') and Gf[k][i] != 0 ):
            if( i == len(Gf[0])-1 ):
                return P.append(i)
            elif( i not in P ):
                P.append(i)
                if(findAugmPathDepth(Gf,i,P) != []):
                    return P
                else: P.pop()
    return []

def findAugmPathBreadth(Gf):
    level = [(0,-1)]
    marked = [(0,-1)]
    C = [0]
    while level != []:
        next_level = []
        for u in level:
            for i in range(len(Gf[0])):
                if( Gf[u[0]][i] != 0 and i not in C ):
                    if( i == len(Gf[0])-1 ):
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

def updateResNet(Gf,mx,P):
    for i in range(np.size(P)-1):
        Gf[P[i]][P[i+1]] -= mx
        Gf[P[i+1]][P[i]] += mx


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
    return np.sqrt( (p[0]-q[0])**2 + (p[1]-q[1])**2 )+1e-8

def histograma(i,P):
    radBins = 12        # Cantidad de cajas radiales (anillos)
    angleBins = 10      # Cantidad de cajas angulares (quesitos)
    firstBin = -2       # Exponente del 2 del primer bin
    PA = np.copy(P)
    """
    X =[0,0]
    for i in range(len(PA)):
        X[0] += 1.0*PA[i][0]/len(PA)
        X[1] += 1.0*PA[i][1]/len(PA)
    """
    # Normaliza la forma con centro en X
    X = P[i]
    Max = max([norm(P[i],X) for i in range(len(P))])
    for i in range(len(PA)):
        PA[i][0] = 1.0*(P[i][0]-X[0])/(Max)
        PA[i][1] = 1.0*(P[i][1]-X[1])/(Max)

    # Inicializa el histograma
    H = np.zeros([radBins,angleBins])

    # Llena los bins (matriz de presencia)
    for j in range(len(PA)):
        if(j != i):
            theta = np.arctan2(PA[j][1],PA[j][0])
            r = norm(PA[i],PA[j])
            if int(np.log2(r)+firstBin) < radBins+ firstBin:
                H[int(np.log2(r)+firstBin)-firstBin][int(1.0*theta/(2*np.pi))] += 1.0
    return H/float(sum(sum(H))

def costFunction(ptInFirstShape,ptInSecondShape, shape1,shape2):
    H1 = histograma(ptInFirstShape,shape1)
    H2 = histograma(ptInSecondShape,shape2)

    C = 0
    for k1 in range(H1):
        for k2 in range(H1[0]):
            C += ((H1[k1][k2]-H2[k1][k2])**2)/(H1[k1][k2]+H2[k1][k2])
    return C/2.0

def normalize(Y):
    X = np.zeros((len(Y),2))
    Q = [0,0]
    for i in range(len(Y)):
        Q[0] += 1.0*Y[i][0]/len(Y)
        Q[1] += 1.0*Y[i][1]/len(Y)
    N = max([norm(Q,Y[i]) for i in range(len(X))])
    for i in range(len(X)):
        X[i][0] = float(Y[i][0]-Q[0])/N
        X[i][1] = float(Y[i][1]-Q[1])/N
    return X

def shapeContext(X, Y):
    C = np.zeros([size,size])

    for i in range(size):
        for j in range(size):
            if i != j:
                C[i][j] = costFunction(i,j,X,Y)
    return C

def shapenorm(X,Y):
    sizeX = len(X)
    sizeY = len(Y)
    size = 0
    if sizeX > sizeY:
        for _ in range(sizeX-sizeY):
            Y = np.append(Y,[2**12,2**12])
        size = sizeX
    elif sizeY > sizeX:
        for _ in range(sizeY-sizeX):
            X = np.append(X,[2**12,2**12])
        size = sizeY
    else: size = sizeX

    C = shapeContext(X,Y)
    A = hungarianAlgorithm(C)
    norm = 0
    for (i,j) in A:
        norm += C[i-1][j-1]
    return [A,norm]







#########################
# Funciones auxililares #
#########################

#Funcion que determina el tipo de vuelta de las rectas OX y OY
def CWorCCW(X = (0,0), O  = (0,0), Y  = (0,0)):
	#Ve si se repiten puntos
	if X == O or X == Y or O == Y: return "Nil"
	#Revisa la tricotomia entre las pendientes de ambas rectas
	if (X[1] - O[1])*(Y[0] - O[0]) < (Y[1] - O[1])*(X[0] - O[0]): return -1
	elif (X[1] - O[1])*(Y[0] - O[0]) > (Y[1] - O[1])*(X[0] - O[0]): return 1
	else: return 0

#Funcion que cambia 2 variables
def swap(a,b):
	#Cambia sin usar una tercera variable
	try:
		a = a+b
		b = a-b
		a = a-b
		return a,b
	#En caso de un overflow en la primera suma
	except:
		aux = a
		a = b
		b = aux
		return a,b

#Funcion que ordena los puntos con bubblesort
def sortPoints(X,Y):
	l = 1
	while l != 0:
		l = 0
		for i in range(len(X)-1):
			if (X[i] > X[i+1]) or (X[i] == X[i+1] and Y[i] > Y[i+1]):
				(X[i],X[i+1]) = swap(X[i],X[i+1])
				(Y[i],Y[i+1]) = swap(Y[i],Y[i+1])
				l += 1

#Funcion que encuentra la distancia de un punto a una linea, el signo determina de que lado
def distPoint_Line(P0,P1,P2, sign = 0):
	if sign == 0:
		return abs( (P2[1]-P1[1])*P0[0] - (P2[0]-P1[0])*P0[1] + P2[0]*P1[1] - P2[1]*P1[0] ) / ((P2[1]-P1[1])**2 + (P2[0]-P1[0])**2)**0.5
	#MENOS: El punto esta CCW la linea, o la dist. sera negativa
	elif sign == -1:
		return ( (P2[1]-P1[1])*P0[0] - (P2[0]-P1[0])*P0[1] + P2[0]*P1[1] - P2[1]*P1[0] ) / ((P2[1]-P1[1])**2 + (P2[0]-P1[0])**2)**0.5
	#MAS: El punto esta CW la linea, o la dist. sera negativa
	elif sign == 1:
		return -( (P2[1]-P1[1])*P0[0] - (P2[0]-P1[0])*P0[1] + P2[0]*P1[1] - P2[1]*P1[0] ) / ((P2[1]-P1[1])**2 + (P2[0]-P1[0])**2)**0.5

#Funcion que encuentra el punto a mayor distancia de una recta
def maxdistP_L(X,Y,A,P1,P2, sign = 0):
	aux = A[0]
	if sign == -1:
		for i in A:
			if distPoint_Line( (X[aux],Y[aux]), P1,P2, sign ) < distPoint_Line( (X[i],Y[i]), P1,P2, sign ):
				aux = i
	else:
		for i in A:
			if distPoint_Line( (X[aux],Y[aux]), P1,P2, sign ) > distPoint_Line( (X[i],Y[i]), P1,P2, sign ):
				aux = i
	return aux

#Funcion que encuentra el angulo entre 3 puntos

def findAngle(P0,P1,P2):
	I = (float('inf'),-float('inf'))
	#if P0[0] in I or P0[1] in I or P1[0] in I or P1[1] in I or P2[0] in I or P2[1] in I:
	#	return 0
	a = (P1[0]-P0[0])**2 + (P1[1]-P0[1])**2
	b = (P1[0]-P2[0])**2 + (P1[1]-P2[1])**2
	c = (P2[0]-P0[0])**2 + (P2[1]-P0[1])**2
	return np.arccos( (a+b-c)/(4*a*b+1e-7)**0.5 )



def grahamScan(X,Y):
	#Ordena los puntos en orden lexicografico
	sortPoints(X,Y)
	#Inicia las listas superior e inferior, resp.
	(Lu,Ll) = ([],[])
	#Mete los primeros 2 puntos en la lista superior
	Lu.append(0)
	Lu.append(1)
	#Este ciclo recorre los puntos que faltan
	for i in range(2,len(X)):
		#En cada iteracion prueba cada punto
		Lu.append(i)
		#Este ciclo revisa si se ha dado una vuelta mala
		while len(Lu) > 2 and CWorCCW( (X[Lu[-3]],Y[Lu[-3]]),(X[Lu[-2]],Y[Lu[-2]]), (X[Lu[-1]],Y[Lu[-1]]) ) == 1:
				Lu.pop(-2)
	#Esto funciona, introduce los ultimos 2 puntos a la lista inferior
	Ll.append(len(X)-1)
	Ll.append(len(X)-2)
	#Repite el proceso anterior, solo que recorriendo en orden contrario
	for i in range(len(X)-3,-1,-1):
		Ll.append(i)
		while len(Ll) > 2 and CWorCCW( (X[Ll[-3]],Y[Ll[-3]]),(X[Ll[-2]],Y[Ll[-2]]), (X[Ll[-1]],Y[Ll[-1]]) ) == 1:
			Ll.pop(-2)
	#Fusiona y regresa las listas (por fines graficos, incluyen los extremos)
	Lu.pop()
	return Lu+Ll

def Chan(X,Y):
	t = 1
	while True:
		#Busca acotar la cantidad de puntos en el convex hull
		m = min(2**(2**t),len(X))
		#Esto es para fines practicos
		if len(X)/m == 1.0*len(X)/m:
			r = int(len(X)/m)
		else:
			r = int(len(X)/m)+1
		#Inicia la lista de convex hull's
		P,CH, R, Aux = [],[],[], []
		for i in range(r):
			#Conforme agrega alementos los convierte en sub convex hull's
			CH.append([])
			#Esta es una lista de listas auxiliares
			R.append([])
			R[i] = [k for k in range(i*m,min((i+1)*m,len(X)))]
			#Realiza Graham Scan en cada pequenha lista
			Aux = grahamScan(X[R[i]],Y[R[i]] )
			#Esto recupera los indices de sus respectivas listas (no solo 1,2,3,4)
			CH[i] = [R[i][Aux[j]] for j in range(len(Aux)-1)]
			#Evita elementos repetidos de cierre de ciclo
			CH[i].pop()
		#Inicializa la busqueda de angulos
		l = 0
		aux = 0
		for i in range(1,len(X)):
			if Y[i] < Y[l]: l = i
		P.append( l )
		P.insert(0, (X[l]+0.1,Y[l]) )
		for k in range(1,m):
			Q = []
			for i in range(r):
				Q.append([])
				if k == 1: aux = max([(findAngle( P[k-1],(X[P[k]],Y[P[k]]),(X[q],Y[q]) ), q) for q in range(i*m,min((i+1)*m,len(X)))])
				else: aux = max([(findAngle( (X[P[k-1]],Y[P[k-1]]),(X[P[k]],Y[P[k]]),(X[q],Y[q]) ), q) for q in range(i*m,min((i+1)*m,len(X)))])
				Q[i] = aux
				#except: Q[i] = ()
				#print Q[i]
			if k == 1: aux = max([(findAngle( P[k-1],(X[P[k]],Y[P[k]]),(X[Q[i][1]],Y[Q[i][1]]) ), Q[i][1]) for i in range(r)] )
			else: aux = max([(findAngle( (X[P[k-1]],Y[P[k-1]]),(X[P[k]],Y[P[k]]),(X[Q[i][1]],Y[Q[i][1]]) ), Q[i][1]) for i in range(r)] )
			P.append( aux[1] )
			if P[k+1] == P[1]:
				P.pop(0)
				return P
		t += 1


def Quickhull(X,Y):
	if X == []:
		return []
	k1,k2 = 0,0
	for i in range(1,len(X)):
		#Cambia punto objetivo si encuentra uno mas a la izquierda
		if X[i] < X[k1]: k1 = i
		if X[i] > X[k2]: k2 = i
	#Inicia las listas superior e inferior, resp.
	Lu,Ll = [k2,k1],[k1,k2]
	#Separa los puntos segun su posicion respecto la linea en 2 subconj.
	Qu, Ql = [],[]
	for k in range(len(X)):
		if CWorCCW( (X[k1],Y[k1]), (X[k],Y[k]), (X[k2],Y[k2]) ) == -1:
			Qu.append(k)
		else:
			Ql.append(k)

	#Resuelve el convex hull superiormente:
	#Indica si se realizan cambios
	it = 1
	#Mientras se realizen cambios
	while it != 0:
		#Reinicia las variables auxiliar a 0
		it = 0
		aux = 0
		#Para cada segmento en el poligono contruido...
		for i in range(len(Lu)-1):
			#Encuentra el punto mas lejano que este fuera
			aux = maxdistP_L(X,Y,Qu,(X[Lu[i]],Y[Lu[i]]), (X[Lu[i+1]],Y[Lu[i+1]]), 1)
			#Si se encuentra fuera...
			if distPoint_Line( (X[aux],Y[aux]), (X[Lu[i]],Y[Lu[i]]), (X[Lu[i+1]],Y[Lu[i+1]]), 1) < 0:
				#Introduce al poligono
				Lu.insert(i+1,aux)
				#Indica que se realizo cambio
				it = 1
				#Quita los puntos en el Hull
				for k in Qu:
					#Quita los puntos interioresal poligono
					try:
						if CWorCCW((X[Lu[i]],Y[Lu[i]]),(X[k],Y[k]),(X[aux],Y[aux])) + CWorCCW((X[aux],Y[aux]),(X[k],Y[k]),(X[Lu[i+2]],Y[Lu[i+2]])) + CWorCCW((X[Lu[i+2]],Y[Lu[i+2]]),(X[k],Y[k]),(X[Lu[i]],Y[Lu[i]])) == -3:
							Qu.remove(k)
					#Quita los puntos sobre el poligono
					except: Qu.remove(k)
	#Resuelve inferiormente el convex hull

	#Vuelve a iniciar la variable auxiliar que indica si hubo cambios
	it = 1
	#Realiza lo mismo, pero con la parte inferior
	while it != 0:
		it = 0
		aux = 0
		for i in range(len(Ll)-1):
			aux = maxdistP_L(X,Y,Ql,(X[Ll[i]],Y[Ll[i]]), (X[Ll[i+1]],Y[Ll[i+1]]), 1)
			if distPoint_Line( (X[aux],Y[aux]), (X[Ll[i]],Y[Ll[i]]), (X[Ll[i+1]],Y[Ll[i+1]]), 1) < 0:
				Ll.insert(i+1,aux)
				it = 1
				for k in Ql:
					try:
						if CWorCCW((X[Ll[i]],Y[Ll[i]]),(X[k],Y[k]),(X[aux],Y[aux])) + CWorCCW((X[aux],Y[aux]),(X[k],Y[k]),(X[Ll[i+2]],Y[Ll[i+2]])) + CWorCCW((X[Ll[i+2]],Y[Ll[i+2]]),(X[k],Y[k]),(X[Ll[i]],Y[Ll[i]])) == -3:
							Ql.remove(k)
					except: Ql.remove(k)
	#Para no repetir puntos innecesarios
	Lu.pop()
	#Une y regresa la lista con los puntos en el convex hull
	return Lu+Ll

#img = cv2.imread('Documentos/Desarrollo/Python/Geometria\\Computacional/Project/Im\'agenes/')

def checkColor(col,minI = [0,0,0], maxI = [15,15,15]):
	if col[0] >= minI[0] and col[0] <= maxI[0]:
		if col[1] >= minI[1] and col[1] <= maxI[1]:
			if col[2] >= minI[2] and col[2] <= maxI[2]:
				return False
	else: return True

def samplePoints(img,N):
	i = 0
	S = []
	while i < N:
		A = 256*np.random.random(2)
		col = img[int(A[0]),int(A[1])]
		if checkColor(col) == True:
			S.append([int(A[0]),int(A[1])])
			i += 1
	return S

def concatenateV(C1,C2,P):
	C1.pop()
	C2.pop()
	y1 = max([[P[C1[i]][1],i] for i in range(len(C1))])
	y2 = min([[P[C2[i]][1],i] for i in range(len(C2))])
	C = C1[y1[1]:]+C1[:y1[1]]+C2[y2[1]:]+C2[:y2[1]]
	C.insert(0,C[0])
	return C

def concatenateH(C1,C2,P):
	C1.pop()
	C2.pop()
	y1 = min([[P[C1[i]][0],i] for i in range(len(C1))])
	y2 = max([[P[C2[i]][0],i] for i in range(len(C2))])
	C = C2[y2[1]:]+C2[:y2[1]]+C1[y1[1]:]+C1[:y1[1]]
	C.insert(0,C[0])
	return C

def sectImage(X,Y,img):
	P = samplePoints(img,N)
	M = []
	C = []
	for i in range(D[0]):
		C.append([])
		for j in range(D[1]):
			C[i].append([])

	Ck = []
	for i in range(D[0]):
		M.append([])
		for j in range(D[1]):
			M[i].append([])
	for k in range(len(P)):
		M[int((P[k][0]*D[0])/256)][int((P[k][1]*D[1])/256)].append([P[k],k])

	#Calcula los convex Hull de cada uno
	for i in range(D[0]):
		for j in range(D[1]):
			O = np.asarray([M[i][j][k][0] for k in range(len(M[i][j]))])
			C[i][j] = Quickhull(O[:,0],O[:,1])
			for k in range(len(C[i][j])):
				C[i][j][k] = M[i][j][C[i][j][k]][1]

	for i in range(D[0]):
		for j in range(1,D[1]):
			C[i][j] = concatenateV(C[i][j-1],C[i][j],P)
	for k in range(D[0]):
		Ck.append(C[k][-1])
	for k in range(1,D[0]):
		Ck[k] = concatenateH(Ck[k],Ck[k-1],P)



	#Ck = concatenateV(C[0][0],C[0][1],P)
	#Ck = concatenateV(C[1][0],C[1][1],P)
	#Ck = concatenateH(C[1][0],C[0][0],P)

	P = np.asarray(P)
	plt.plot(P[:,0],P[:,1], '.b')
	plt.plot(P[Ck[-1],0],P[Ck[-1],1],'-r')
	plt.show()

	return [[P[Ck[-1][i]][0],P[Ck[-1][i]][1]] for i in range(len(Ck[-1]))]





img1 = cv2.imread('foco1[].jpg')
img2 = cv2.imread('foco2[].jpg')
"""
P = samplePoints(img,N)

P = np.asarray(P)

C = Chan(P[:,0],P[:,1])

plt.plot(P[:,0],P[:,1], '.b')
plt.plot(P[C,0],P[C,1],'-r')
plt.show()
"""
S1 = sectImage(X,Y,img1)
S1 = np.asarray(S1)

S2 = sectImage(X,Y,img2)
S2 = np.asarray(S2)

Res = shapenorm([S1[i] for i in range(0,len(S1),scale)],[S2[i] for i in range(0,len(S2),scale)])

A = Res[0]
dist = Res[1]
Pl = []

for i in range(len(A)):
    Pl.append((S1[scale*A[i][1]-scale][0],S1[scale*A[i][1]-scale][1]))
    Pl.append((S2[scale*A[i][0]-scale][0],S2[scale*A[i][0]-scale][1]))
for i in range(0,len(Pl),2):
    plt.subplot(1,1,1),plt.plot([Pl[i][0],Pl[i+1][0]],[Pl[i][1],Pl[i+1][1]], '--')

plt.subplot(1,1,1),plt.plot(S1[:,0],S1[:,1],'-b')
plt.subplot(1,1,1),plt.plot(S2[:,0],S2[:,1],'-r')

print('La distancia entre las imagenes es: ', dist)
plt.show()
