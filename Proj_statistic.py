#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
1
"""
import numpy as np
from astropy.io import ascii
import networkx as nx
#from networkx.generators.degree_seq import expected_degree_graph
import matplotlib.pyplot as plt
import random
Ns=[5,10,20,30,40,50]

Ds=[]
As=[]
Ss=[]
ans=[]
er1=[]
er2=[]
cc=-1
ef_av1pl=[]
ef_av2pl=[]
for N in (Ns):
    #Ds=[]
    #As=[]
    #Ss=[]
    #ans=[]
    er1=[]
    er2=[]
    for ensemble in range(100):
        print(ensemble,'\n and N=',N)
        #cc=cc+1
        #G=nx.gnm_random_graph(N,N) #code creating G here
        
        G=nx.random_powerlaw_tree(N, gamma=3, seed=None, tries=10000)
        #for (u,v,w) in G.edges(data=True):
        #    w['weight'] = random.randint(0,1)
        
        #for (u, v) in G.edges():
        #    G.edges[u,v]['weight'] = random.randint(0,1)
        #pos=nx.spring_layout(G,iterations=0)
        #nx.draw_networkx(G)
        #labels = nx.get_edge_attributes(G,'weight')
        #nx.draw_networkx_edge_labels(G,edge_labels=labels)
        A = nx.to_numpy_matrix(G)
        A=np.squeeze(np.asarray(A))
        
        #print(A)
        sigma=1
        S=np.zeros((np.shape(A)[0],20))
        
        for i in range (np.shape(S)[0]):
            S[i][0]=random.uniform(0,20)
        
        for i in range (np.shape(S)[0]):
            S[i][1]=S[i][0]+random.uniform(-1, 1)
        
        def product(A,D):
            p = []
            for i in range (np.shape(A)[0]):
                sss = 0
                for j in range(np.shape(A)[0]):
                    sss = sss + A[i][j]*D[j]
                p.append(sss)
            return p
            
        
        
        D2=[]
        for t in range (2,np.shape(S)[1]):
            D=[]
            for i in range (np.shape(S)[0]):
                K=(S[i][t-2])-S[i][t-1]
                D.append(K)
            D2.append(D)
            for i in range (np.shape(S)[0]):
                #product(A,D).reverse()
                S[i][t]=S[i][t-1]+(product(A,D)[i])*sigma#+random.uniform(-0.2, 0.2)
                if (S[i][t])>=20:
                    S[i][t]=20
                if (S[i][t])<=0:
                    S[i][t]=0
        #Ds.append(D2)
        #As.append(A)
        #Ss.append(S)
        D2=np.transpose(D2)
        TH=1/5
        
        def EV(b, n): 
            prb = 1 / n 
            sum = 0
            for i in range(0, n): 
                sum += (b[i] * prb)   
            return float(sum)
        
        #A=As[cc]
        DD=[]
        for i in range(len(D2)):
            y=[]
            for j in range(len(D2[0])):
                if (j+1)<len(D2[0]):
                    x=D2[i][j+1]-D2[i][j]
                    y.append(x)
            DD.append(y)
        analizer=[]
        student=[i for i in range (len(DD))]
        for i in range(len(student)):
           x=[]
           for j in range(len(DD)):
               y=[]
               for k in range(len(DD[0])):
                   z=round(abs(DD[i][k])-abs(DD[j][k]),2)
                   y.append(z)
               x.append(y)
           analizer.append(x)
        matrix=np.zeros((len(D2),len(D2)))  
        for i in range(len(D2)):
            for j in range(len(D2)):
                matrix[i][j]=round(abs(EV(analizer[i][j],len(analizer[i][j]))),2)
        a=np.zeros((len(D2),len(D2)))
        for i in range(len(D2)):
            for j in range(len(D2)):
                if matrix[i][j]<=np.mean(matrix)*TH:
                    a[i][j]=1
                    if (D2[i][1]==0 or D2[j][1]==0)and (D2[i][2]==0 or D2[j][2]==0) :
                        a[i][j]=0
                        a[j][i]=0
                    if i==j or j==i:
                        a[i][j]=0
                        a[j][i]=0
                    else:
                        pass
                else:
                    a[i][j]=0
        
        #ans.append(a)
        list_A=[]
        list_a=[]
        for i in range(len(A[0])):
            p=0
            for j in range(len(A[0])):
                p+=A[i][j]
            list_A.append(p)
        for i in range(len(a[0])):
            p=0
            for j in range(len(a[0])):
                p+=a[i][j]
            list_a.append(p)
        error=[(list_A[i])-(list_a[i]) for i in range(len(list_a))]
        for l in range(len(error)):
            if error[l]<0:
                error[l]=(-1)*error[l]
        #print(error,'\n')
        ef=(np.average(error))/(N-1)
        ef=abs(ef-1)
        er1.append(ef)
        
        aa=np.zeros((len(a),len(a)))
        q=np.zeros((len(A),len(A)))
        for i in range(len(A)):
            for j in range(len(A)):
                q[i,j]=A[i][j]
        for i in range(len(a)):
            for j in range(len(a)):
                aa[i,j]=a[i][j]
        error=np.zeros((len(A[0]),len(A[0])))
        for i in range(len(A[0])):
            for j in range(len(A[0])):     
                error[i,j]=q[i,j]-aa[i,j]
        for i in range(N):
            for j in range (N):
                if error[i][j]<0:
                    error[i][j]=error[i][j]*(-1)
        ef=np.sum(error)/((N**2)-N)
        er2.append(ef)
    ef_av1pl.append(np.average(er1))
    ef_av2pl.append(np.average(er2))
        #print(error) 
# plt.plot(Ns,ef_av1pl)
# plt.xlabel('N')
# plt.ylabel('efficiency')
# plt.title('efficiency metthod 1 in 100 ensemble for each N (PL network)')
# plt.figure()
# plt.plot(Ns,ef_av2pl)
# plt.xlabel('N')
# plt.ylabel('efficiency')
# plt.title('efficiency metthod 2 in 100 ensemble for each N (PL network)')
# plt.figure()
# plt.plot(Ns,ef_av1pl,label='method 1')
# plt.plot(Ns,ef_av2pl,label='method 2')
# plt.xlabel('N')
# plt.ylabel('efficiency')
# plt.title('efficiency metthod 1,2 in 100 ensemble for each N (PL network)')
# plt.legend()

# plt.plot(Ns,ef_av1pl,label='method 1 PL')
# plt.plot(Ns,ef_av1,label='method 1 ER')
# plt.plot(Ns,ef_av2pl,label='method 2 PL')
# plt.plot(Ns,ef_av2,label='method 2 ER')
# plt.xlabel('N')
# plt.ylabel('efficiency')
# plt.title('efficiency metthod 2 in 100 ensemble for each N (PL,ER network)')
# plt.legend()




