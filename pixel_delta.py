from math import *
import matplotlib.pyplot as plt
import numpy as np
	
# alpha is half of vertical fov
# sigma is angle between laser plane and camera center vector
# phi is angle between camera center vector and target pixel vector
# d is distance between laser plane and camera
# h is vertical pixel count
alpha = pi/6
sigma = pi/3
d = 5
h = 1000

def get_phi(y, h, alpha):
	# get phi of pixel row y
	return atan(tan(alpha)*(2*y/h - 1))

def crunch(alpha, sigma, phi, d, h):
	# depth delta approximation
	da = 2*d*tan(alpha)*cos(phi)/(h*sin(sigma+phi)**2)
	# te is tan(epsilon)
	# epsilon is angle increment in phi across one pixel
	te = (2*tan(alpha)*cos(phi)**2/(h + 2*tan(alpha)*sin(phi)*cos(phi)))
	# exact depth delta
	#de = d/sin(sigma+phi)*te/(sin(sigma+phi) + te*sin(sigma+phi))
	de = 2*d*tan(alpha)*cos(phi)*(cos(phi) - sin(phi)*te)/(h*sin(sigma+phi)*(sin(sigma+phi) + cos(sigma+phi)*te))
	return {'da':da, 'te':te, 'de':de, 'ea':de-da, 'er':(de-da)/de}

def nprint(phi_text, phi):
	nums = crunch(alpha, sigma, phi, d, h)
	print(phi_text)
	print(f"Delta approx.:    {nums['da']}")
	print(f"Tan(eps):         {nums['te']}")
	print(f"Delta exact:      {nums['de']}")
	print(f"Error (absolute): {nums['ea']}")
	print(f"Error (relative): {nums['er']}")
	print("")

nprint('Maximal phi', alpha)
nprint('Zero phi', 0)
nprint('Minimal phi', -alpha)

print(f"phi at y = 50 is {get_phi(50, h, alpha)}")

# make data
xpoints = []
dapoints = []
depoints = []
epoints = []
aspan = 0
espan = 0
tspan = d*(1/tan(sigma-alpha) - 1/tan(sigma+alpha))
for y in range(h):
	xpoints.append(y)
	cronch = crunch(alpha, sigma, get_phi(y, h, alpha), d, h)
	dapoints.append(cronch['da'])
	aspan += cronch['da']
	depoints.append(cronch['de'])
	espan += cronch['de']
	epoints.append(cronch['er'])

print('approximate segment sum', aspan)
print('exact segment sum', espan)
print('expected sum', tspan)

xpoints = np.array(xpoints)
dapoints = np.array(dapoints)
depoints = np.array(depoints)
plt.plot(xpoints, dapoints, label = 'approx')
plt.plot(xpoints, depoints, label = 'exact')
plt.legend()
plt.show()
plt.plot(xpoints, epoints, label = 'error')
plt.legend()
plt.show()
