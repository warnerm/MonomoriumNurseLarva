import numpy as np
import matplotlib.pyplot as plt
import math as m

Genes = np.loadtxt("/Users/michaelwarner/MockGenes.txt",skiprows=1,dtype='string')
Conns = np.loadtxt("/Users/michaelwarner/MockConnections.txt",skiprows=1,dtype='int')

plt.figure(figsize=(6,6),facecolor='white')
plt.subplot(111)

n=3
circSize = 2000
textOffset = 0.05
r = 0.35
center = [[0.5,0,1],[0.1,1,1]]
colors = ["#000000", "#E69F00", "#56B4E9"]


def Addcircles():
	for i in range(3):
		AddOneType(center[0][i],center[1][i],colors[i])
	return

def AddOneType(cx, cy,color):
	P = CircleCoord (cx,cy)
	S = np.linspace(circSize,circSize,n*2)
	scat = plt.scatter(P[:,0],P[:,1],s=S,lw = 0.5,marker="^",color=color)
	return

def CircleCoord ( cx, cy):
	X1 = np.fromfunction(lambda i, j: cx - r*np.cos((i*np.pi/(n))), (n,1),dtype=float)
	X2 = np.fromfunction(lambda i, j: cx - r*np.cos(((n-i)*np.pi/(n))), (n,1),dtype=float)
	X = np.concatenate((X1,X2),axis=0)
	Y1 = np.sqrt(r ** 2 - (X1 - cx) ** 2) + cy
	Y2 = -np.sqrt(r ** 2 - (X1 - cx) ** 2) + cy
	Y = np.concatenate((Y1,Y2),axis=0)
	P = np.hstack((X,Y))
	return P

def Addnames ():
	for i in range(3):
		addLabels(i);
	return

def addLabels (Type):
	P = CircleCoord (center[0][Type],center[1][Type])
	for i in range(6):
		plt.text(P[i][0],P[i][1]-textOffset,Genes[i][Type],
			fontsize=16,horizontalalignment='center',
			verticalalignment='center',color="w")
	return

Addcircles(); #Add circles for genes
Addnames(); #Add gene names
#Addlinks(); #Add arrows
plt.axis('off')
plt.show()





