import math

def exact(x, n, h, l) :
	d = 1.0/n # dx
	a = d*d # Area of crosssection
	p = d*4 # boundary length
	m = math.sqrt( (h*p)/(l*a) )
	c = math.cosh(m*(1-x))/math.cosh(m)
	print x, c

# main
n = 5    # n is num. of division
h = 12    # h is coef. of heat transfer
l = 50    # l is coef. of heat conduction
d = 1.0/n # d is dx
st = d*0.5

for i in range(n) :
	exact(st+i*d, n, h, l)