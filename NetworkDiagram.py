import numpy as np
import matplotlib.pyplot as plt
import math as m

plt.figure(figsize=(6,6),facecolor='white')
plt.subplot(111)

n=100
cx = cy = 0.5
r = 0.25
center = [[0.5,0.33,0.67],[0.33,0.67,0.67]]



def Addcircles():
	P = CircleCoord (cx,cy,r,n)
	S = np.linspace(50,50,n*2)
	scat = plt.scatter(P[:,0],P[:,1],s=S,lw = 0.5)
	plt.axis('off')
	plt.show()
	return



def CircleCoord ( cx, cy, r, n):
	X1 = np.fromfunction(lambda i, j: cx - r*np.cos((i*np.pi/(n))), (n,1),dtype=float)
	X2 = np.fromfunction(lambda i, j: cx - r*np.cos(((n-i)*np.pi/(n))), (n,1),dtype=float)
	X = np.concatenate((X1,X2),axis=0)
	Y1 = np.sqrt(r ** 2 - (X1 - cx) ** 2) + cy
	Y2 = -np.sqrt(r ** 2 - (X1 - cx) ** 2) + cy
	Y = np.concatenate((Y1,Y2),axis=0)
	P = np.hstack((X,Y))
	return P


Addcircles(); #Add circles for genes
#Addnames(); #Add gene names
#Addlinks(); #Add arrows






