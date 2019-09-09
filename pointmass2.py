from matplotlib import pyplot as plt
from math import sqrt

m1 = 0.1
m2 = 0.1
dis = 0.1
tol = 1.0
def alphax(kx, ky):	
	return kx - (m1*(kx-dis)/((kx-dis)**2 + (ky)**2)) - (m2*(kx+dis)/((kx+dis)**2 + (ky)**2))

def alphay(kx, ky):	
	return ky - (m1*ky/((kx-dis)**2 + (ky)**2)) - (m2*ky/((kx+dis)**2 + (ky)**2))

def cross(a, b):
	return a[0]*b[1] - a[1]*b[0]

def dist(a,b):
	return sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)

def contains(aix, bix, cix, x):	
	da = (x[0]-aix[0], x[1]-aix[1])
	db = (x[0]-bix[0], x[1]-bix[1])
	dc = (x[0]-cix[0], x[1]-cix[1])

	if min(dist(da,db), dist(db,dc), dist(dc,da)) > tol:
		return False


	if (cross(da, db)>=0 and cross(db, dc)>=0 and cross(dc, da)>=0) or (cross(da, db)<=0 and cross(db, dc)<=0 and cross(dc, da)<=0):
		return True
	else:
		return False

points = []	
def findimagesfor(bx,by):	
	xstart = -4.0
	ystart = -4.0
	xstep = -xstart*2.0/1024
	ystep = -ystart*2.0/1024
	xstart = -3.98
	ystart = -3.98
	kxvalues = [xstart+i*xstep for i in range(1024)]
	kyvalues = [ystart+i*ystep for i in range(1024)]

	for kx in kxvalues:
		for ky in kyvalues:
			a = (kx,ky)
			b = (kx+xstep, ky)
			c = (kx+xstep, ky+ystep)
			d = (kx, ky+ystep)
			aix, aiy = (alphax(a[0],a[1]), alphay(a[0],a[1]))
			bix, biy = (alphax(b[0],b[1]), alphay(b[0],b[1]))
			cix, ciy = (alphax(c[0],c[1]), alphay(c[0],c[1]))
			dix, diy = (alphax(d[0],d[1]), alphay(d[0],d[1]))

			if(contains((aix,aiy), (bix,biy), (cix, ciy), (bx,by))):
				points.append((a[0] , a[1]))				

			elif(contains((aix,aiy), (cix,ciy), (dix, diy), (bx,by))):
				points.append((a[0] , a[1]))				

	print
	print(points)
	x = [p[0] for p in points]
	y = [p[1] for p in points]	
	
	print
	plt.scatter(x, y, color='blue')
	plt.scatter([bx,],[by,], color='red')
	plt.scatter([dis, -dis], [0, 0], color='green')	
	plt.show()

	return points

findimagesfor(0.0, 0.5)