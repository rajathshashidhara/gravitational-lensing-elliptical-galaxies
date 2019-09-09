from matplotlib import pyplot as plt
from math import sqrt

#tol = 0.05
def alphax(kx, ky):
	m = 1
	return kx - (m*kx/((kx)**2 + (ky)**2))

def alphay(kx, ky):
	m = 1
	return ky - (m*ky/((kx)**2 + (ky)**2))

def cross(a, b):
	return a[0]*b[1] - a[1]*b[0]

def contains(aix, bix, cix, x):
	da = (x[0]-aix[0], x[1]-aix[1])
	db = (x[0]-bix[0], x[1]-bix[1])
	dc = (x[0]-cix[0], x[1]-cix[1])

	if (cross(da, db)>=0 and cross(db, dc)>=0 and cross(dc, da)>=0) or (cross(da, db)<=0 and cross(db, dc)<=0 and cross(dc, da)<=0):
		return True
	else:
		return False

def dis(a,b):
	return sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)

def findimagesfor(bx,by):
	points = []

	xstart = -2.0
	ystart = -2.0
	xstep = -xstart*2.0/200
	ystep = -ystart*2.0/200
	xstart = -2.00
	ystart = -2.00
	kxvalues = [xstart+i*xstep for i in range(200)]
	kyvalues = [ystart+i*ystep for i in range(200)]

	for kx in kxvalues:
		for ky in kyvalues:
			a = (kx,ky)
			b = (kx+xstep, ky)
			c = (kx+xstep, ky+ystep)
			d = (kx, ky+ystep)
			try:
				aix, aiy = (alphax(a[0],a[1]), alphay(a[0],a[1]))
				bix, biy = (alphax(b[0],b[1]), alphay(b[0],b[1]))
				cix, ciy = (alphax(c[0],c[1]), alphay(c[0],c[1]))
				dix, diy = (alphax(d[0],d[1]), alphay(d[0],d[1]))
			except:
				continue

			if(contains((aix,aiy), (bix,biy), (cix, ciy), (bx,by))):
				points.append(((a[0])/1.0 , (a[1])/1.0   ))
			elif(contains((aix,aiy), (cix,ciy), (dix, diy), (bx,by))):
				points.append(( (a[0])/1.0 , (a[1])/1.0   ))			


	print(points)
	x = [i[0] for i in points]
	y = [i[1] for i in points]
	plt.scatter(x,y, color='blue')
	plt.scatter([0],[0], color='green')
	plt.scatter([bx,],[by,], color='red')
	plt.show()

findimagesfor(0.0, 0.0)
print
