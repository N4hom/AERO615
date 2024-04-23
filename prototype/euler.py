import numpy as np
import matplotlib.pyplot as plt

np.set_printoptions(precision=5, suppress=True)

def delta2(q , i ):
	
	return  q[i + 1] - 2 * q[i] + q[i - 1]

def interpolateLeft(q , i ):
	
	return 0.5 * (q[i] + q[i - 1])

def interpolateRight(q, i ):
	return 0.5 * (q[i] + q[i + 1])

class variable(object):
	"""docstring for field"""
	def __init__(self, name : str  , Nx, Ny , top, bottom, left, right ):
		
		self.name = name

		self.Nx = Nx
		self.Ny = Ny

		
		# Define coefficient matrix for the variables to be solved
		self.a = coeff(np.zeros((Nx,Ny)) , np.zeros((Nx,Ny)) , np.zeros((Nx,Ny)) , np.zeros((Nx,Ny)), np.zeros((Nx,Ny)), np.zeros((Nx,Ny)))
		
		self.icMax = Nx + 3
		self.jcMax = Ny + 3
		self.field = np.zeros((Nx+4,Ny+4))
		self.field[2:-2] = value


		self.top = top
		self.bottom = bottom
		self.left = left
		self.right = right
		

		# Place-holders for the line-by-line procedure
		# self.Axs = [None]*Ny
		# self.Ays = [None]*Nx
		# self.bxs = [None]*Ny
		# self.bys = [None]*Nx

class variable1D(object):
	"""docstring for field"""
	def __init__(self, name : str  , Nx , left, right , value = 0.0):
		
		self.name = name

		self.Nx = Nx
		
		self.D = coeff(None , None , np.zeros(Nx) , np.zeros(Nx), np.zeros(Nx), None)

		
		# Define coefficient matrix for the variables to be solved
		self.R = coeff(n = None, s = None , e = np.zeros(Nx) , w = np.zeros(Nx), p = np.zeros(Nx), b = None)
		
		self.iMax = Nx + 3
		# self.jMax = Ny + 1
		self.field = np.zeros((Nx+4))
		self.field[2:-2] = value
		self.flux_f = np.zeros(Nx+4)
		self.flux_g = np.zeros(Nx+4)

		self.residual = []

		
		self.left = left
		self.right = right
		
	def __getitem__(self, i):
		return self.field[i]

	def __setitem__(self, i, val):
		self.field[i] = val

	def print(self):
		print(self.name , np.array_str(self.field, max_line_width=np.inf, precision=5, suppress_small=True))
		return 

	def printFlux(self):
		print(self.name + " flux: " , np.array_str(self.flux_f, max_line_width=np.inf, precision=5, suppress_small=True))
		return 

	

	def computeResidualij(self , i):
		ic = i + 2
		self.R.w[i] =   interpolateLeft(self.flux_f , ic) * 0.1  # face length is hardcoded
		self.R.e[i] =   interpolateRight(self.flux_f , ic) * 0.1

		# if self.name == "rho":
		# 	print("interpolateLeft ",  interpolateLeft(self.flux_f , ic) )
		# 	print("interpolateRight ", interpolateRight(self.flux_f , ic) )


		self.R.p[i] = self.R.e[i] - self.R.w[i]

	def computeResidual(self):
		
		# if self.name == "rho":
		# 	self.printFlux()

		for i in range(len(self.R.p)):
			self.computeResidualij(i)


	def computeDissipation(self , lambdaFace , s2face , s4face):
		
		deltaY = 0.1
		# print("lambda west " , lambdaFace.w)

		for i in range(self.Nx):
			deltaqw  = self.field[i ] - self.field[i - 1]
			delta3qw = delta2(self.field , i) - delta2(self.field , i - 1)
			deltaqe  = self.field[i + 1] - self.field[i ]
			delta3qe = delta2(self.field , i + 1) - delta2(self.field , i )
			self.D.w[i] = s2face.w[i] * deltaY * lambdaFace.w[i] * deltaqw - s4face.w[i] * deltaY * lambdaFace.w[i] * delta3qw
			self.D.e[i] = s2face.e[i] * deltaY * lambdaFace.e[i] * deltaqe - s4face.e[i] * deltaY * lambdaFace.e[i] * delta3qw
			self.D.p[i] = self.D.e[i] - self.D.w[i]

		# print(self.D.w)
		# print(self.D.e)
		# print(self.D.p)


	def computeFluxij(self , u , p , i):
		
		if self.name == "rhoU":
			self.flux_f[i] = self.field[i] * u[i] + p[i]
		elif self.name == "rhoE":
			self.flux_f[i] = self.field[i] * u[i] + p[i] * u[i]
		else :
			self.flux_f[i] = self.field[i] * u[i]
	
	
	def computeFlux(self , u, p):
		
		if self.name == "rhoU":
			self.flux_f = self.field * u + p
		elif self.name == "rhoE":
			self.flux_f = self.field * u + p * u
		else :
			self.flux_f = self.field * u



class coeff(object):

	#This class is a place holder that improves readability.
	def __init__(self, n = None, s = None, w = None, e = None , p =  None , b = None):
		self.n = n
		self.s = s
		self.w = w
		self.e = e

		self.p = p

		self.b = b


# Nx and Ny are the number of control volumes along x and y for the pressure. 

class Euler(object):
	"""2D diffusion class"""
	def __init__(self,  Nx, Ny, Lx , Ly , Mach , pInf , pRatio , rhoInf , gamma = 1.4 ):
		self.Nx = Nx
		self.Ny = Ny
		self.Lx = Lx
		self.Ly = Ly
		self.Deltax = self.Lx/self.Nx  # Deltax and Deltay are default sizes of the CVs
		self.Deltay = self.Ly/self.Ny


		# maximum index with ghost cells
		self.icMax = self.Nx + 3

		self.alpha1 = 0.25
		self.alpha2 = 0.333
		self.alpha3 = 0.75
		self.alpha4 = 1.

		self.pInf = pInf
		self.rhoInf = rhoInf
		self.pRatio = pRatio
		self.M = Mach
		self.cInf = np.sqrt(1.4 * self.pInf / self.rhoInf)
		self.uInf = self.M * self.cInf
		self.gamma = gamma
		self.gamma_1 = gamma - 1

		self.s2 = np.zeros(Nx+4)
		self.s2face = coeff( n = None , s= None , w = np.zeros(self.Nx ) , e = np.zeros(self.Nx))
		self.s4face = coeff( n = None , s= None , w = np.zeros(self.Nx) , e = np.zeros(self.Nx))
		self.lambdaFace = coeff(n = None , s = None , e = np.zeros(Nx + 4) , w = np.zeros(Nx + 4))



		self.rho = variable1D  ("rho"   , Nx , 0.0 , 0.0 , self.rhoInf * 0.1 )
		self.rhoU = variable1D ("rhoU"  , Nx , 0.  , 0.  , self.rhoInf * self.uInf  )  
	#	self.rhoV = variable1D ("rhoV"  , Nx , 1.  , 0.  ,  )
		self.p = variable1D    ("p"     , Nx , 0.  , 0.  , self.pInf * self.pRatio) 
		self.rhoE = variable1D ("rhoE"  , Nx , 0.  , 0.  , self.pInf * self.pRatio/self.gamma_1 + 0.5 * self.rhoInf * self.uInf * self.uInf )
		self.u = variable1D    ("u"     , Nx , 0.  , 0.  , self.uInf )
		self.c = variable1D    ("c"     , Nx , 0.  , 0.  , self.cInf ) # assuming p and rho are inf
		self.M = variable1D    ("M"     , Nx , 0.  , 0.  , self.M)
		self.Riem1 = variable1D("Riem1" , Nx , 0.  , 0.  , self.uInf + 2/self.gamma_1 * self.cInf )
		self.Riem2 = variable1D("Riem2" , Nx , 0.  , 0.  ,  )

	
	def correctInlet(self):
		Riem1 = self.Riem1
		Riem2 = self.Riem2
		p     = self.p
		rho   = self.rho
		rhoU  = self.rhoU
		rhoE  = self.rhoE
		c     = self.c
		u     = self.u
		gamma_1 = self.gamma_1
		gamma = self.gamma
		uInf = self.uInf
		cInf = self.cInf
		M = self.M
		pInf = self.pInf

		print("Correct inlet ")

		Riem1[1] = uInf + 2/gamma_1 * cInf
		Riem2[1] = u[2] - 2/gamma_1 * c[2]  # Riem2[1] = Riem2[2]  

		u[1] = 0.5 * (Riem1[1] + Riem2[1])
		u[0] = u[1]
		
		c[1] = 0.25 * gamma_1 * (Riem1[1] - Riem2[1])
		c[0] = c[1]

		M[1] = u[1]/c[1]
		M[0] = M[1]

		p[1] = pInf/(1 + 0.5 * gamma_1 * M[1]**2)**(gamma/gamma_1)
		p[0] = p[1]

		rho[1] = gamma * p[1]/c[1]**2
		rho[0] = rho[1]

		rhoU[1] = rho[1] * u[1]
		rhoU[0] = rhoU[1]

		rhoE[1] = p[1]/gamma_1 + 0.5 * rho[1] * u[1]**2
		rhoE[0] = rhoE[1]


	def correctOutlet(self):
		Riem1 = self.Riem1
		Riem2 = self.Riem2
		p     = self.p
		rho   = self.rho
		rhoU  = self.rhoU
		rhoE  = self.rhoE
		c     = self.c
		u     = self.u
		gamma_1 = self.gamma_1
		gamma = self.gamma
		uInf = self.uInf
		cInf = self.cInf
		M = self.M
		pInf = self.pInf
		pRatio = self.pRatio
		icMax = self.icMax

		print("Correct outlet")
		
		rho[icMax - 1] = 2*rho[icMax - 2] - rho[icMax - 3]
		rho[icMax    ] = rho[icMax - 1]


		rhoU[icMax - 1] = 2*rhoU[icMax - 2] - rhoU[icMax - 3]
		rhoU[icMax    ] = rhoU[icMax - 1]

		u[icMax - 1] = rhoU[icMax - 1]/rho[icMax - 1]
		u[icMax    ] = u[icMax - 1]

		rhoE[icMax - 1] = pInf*pRatio/gamma_1 + 0.5 * rho[icMax - 1] * u[icMax - 1]**2
		rhoE[icMax    ] = rhoE[icMax - 1]

		p[icMax - 1] = pInf*pRatio
		p[icMax    ] = pInf*pRatio

		c[icMax - 1] = np.sqrt(gamma * p[icMax - 1] / rho[icMax - 1])
		c[icMax    ] = c[icMax - 1]



	def correctFields(self):
		Riem1 = self.Riem1
		Riem2 = self.Riem2
		p     = self.p
		rho   = self.rho
		rhoU  = self.rhoU
		rhoE  = self.rhoE
		c     = self.c
		u     = self.u
		gamma_1 = self.gamma_1
		gamma = self.gamma
		uInf = self.uInf
		cInf = self.cInf
		M = self.M
		pInf = self.pInf

		# pressure and speed of sound shoul be updated only at the interior points
		u.field[2:-2] = rhoU.field[2:-2]/rho.field[2:-2]
		p.field[2:-2] = gamma_1 * (rhoE.field[2:-2] - 0.5 * rho.field[2:-2] * u.field[2:-2] ** 2)
		c.field[2:-2] = np.sqrt(gamma * p.field[2:-2] / rho.field[2:-2])


	def computeSwitches(self):
		
		print("computeSwitches ")

		nu2 = 0
		nu4 = 0.001
		p = self.p

		for i in range(1,len(p.field)-1):
			self.s2[i] = 0.1 * abs(delta2(p , i))/(p[i + 1] + 2*p[i] + p[i - 1]) 

		print("s2 " ,self.s2)
		# Interpolate switches

		for i in range(self.Nx):

			ic = i + 2
			
			self.s2face.w[i] = interpolateLeft(self.s2 , ic) 
			self.s2face.e[i] = interpolateRight(self.s2 , ic) 

			self.s4face.w[i] = max(0 , nu4 - self.s2face.w[i])
			self.s4face.e[i] = max(0 , nu4 - self.s2face.e[i])

			
	def computeEigen(self ):
		print("computeEigen")
		

		for i in range(self.Nx):
			ic = i + 2
			self.lambdaFace.w[ic] = interpolateLeft(self.u.field , ic) + interpolateLeft(self.c.field ,ic)
			self.lambdaFace.e[ic] = interpolateRight(self.u.field , ic) + interpolateRight(self.c.field , ic)



		# print("lambda west " , self.lambdaFace.w)

		# for i in range(self.Nx):
		# 	deltaqw  = field[i ] - field[i - 1]
		# 	delta3qw = delta2(field , i) - delta2(field , i - 1)
		# 	deltaqe  = field[i + 1] - field[i ]
		# 	delta3qe = delta2(field , i + 1) - delta2(field , i )
		# 	field.D.w[i] = self.s2face.w[i] * deltaY * self.lambdaFace.w[i] * deltaqw - self.s4face.w[i] * deltaY * self.lambdaFace.w[i] * delta3qw
		# 	field.D.e[i] = self.s2face.e[i] * deltaY * self.lambdaFace.e[i] * deltaqe - self.s4face.e[i] * deltaY * self.lambdaFace.e[i] * delta3qw
		# 	field.D.p[i] = field.D.e[i] - field.D.w[i]

		# print(field.D.w)
		# print(field.D.e)
		# print(field.D.p)

	def computeFluxes(self):
		
		print("Compute fluxes")
		self.rho.computeFlux(self.u.field , None)
		self.rhoU.computeFlux(self.u.field , self.p.field)
		self.rhoE.computeFlux(self.u.field , self.p.field)
		

	def RungeKutta(self , variable , i):
		
		alpha1 = self.alpha1
		alpha2 = self.alpha2
		alpha3 = self.alpha3
		alpha4 = self.alpha4
		
		dt = 1e-6

		ic = i + 2
		variable0 = variable[ic].copy()

		# if variable.name == "rho":
		# 	print("variable.R.p[i] " , variable.R.p[i])
		# 	print("variable.D.p[i] " , variable.D.p[i])
		# 	print("variable[ic] before" , variable.field[ic])
		variable[i] = variable0 - alpha1 * dt / (self.Lx/self.Nx) * (variable.R.p[i]  - variable.D.p[i])
		variable.computeFluxij(self.u , self.p , ic)
		variable.computeResidualij(i)
		
		variable[i] = variable0 - alpha2 * dt / (self.Lx/self.Nx) * (variable.R.p[i]  - variable.D.p[i])
		variable.computeFluxij(self.u , self.p , ic)
		variable.computeResidualij(i)
		
		variable[i] = variable0 - alpha3 * dt / (self.Lx/self.Nx) * (variable.R.p[i]  - variable.D.p[i])
		variable.computeFluxij(self.u , self.p , ic)
		variable.computeResidualij(i)
		
		variable[ic] = variable0 - alpha4 * dt / (self.Lx/self.Nx) * (variable.R.p[i]  - variable.D.p[i]) # Only one step for testing

		# if variable.name == "rho":
		# 	print("variable[ic] after" , variable[ic])
		

		# variable[i] = variable0 - alpha3 * dt / (self.Lx/self.Nx) * variable.R.p[i]
		# variable[i] = variable0 - alpha4 * dt / (self.Lx/self.Nx) * variable.R.p[i]

		


	def solve(self):

		self.correctInlet()
		self.correctOutlet()
		self.computeFluxes()

		

		self.computeSwitches()
		self.computeEigen()
		
		

		it = 0
		N  = 5000
		while it < N:

			print("-----------------------------")
			print("Iteration " , it)

			self.rho.print()
			self.rhoU.print()
			self.rhoE.print()
			self.p.print()

			print()
			self.rho.printFlux()
			self.rhoU.printFlux()
			self.rhoE.printFlux()
			
			self.computeSwitches()
			print()
			self.computeEigen()
			print()
			
			self.rho.computeDissipation(self.lambdaFace , self.s2face , self.s4face)
			self.rhoU.computeDissipation(self.lambdaFace , self.s2face , self.s4face)
			self.rhoE.computeDissipation(self.lambdaFace , self.s2face , self.s4face)
			
			print()
			
			self.rho.computeResidual()
			self.rhoU.computeResidual()
			self.rhoE.computeResidual()
		
			# self.rho.print()
			# self.rhoU.print()
			# self.rhoE.print()

			print()

			self.correctFields()
			self.computeFluxes()
			
			
			print()

			#Correct boundary
			self.correctInlet()
			self.correctOutlet()
			
			# self.rho.print()
			# self.rhoU.print()
			# self.rhoE.print()

			print()

			rho0 = self.rho.field.copy()

			# loop over internal field
			for i in range(len(self.rho.R.p)):
				# print(i)
				# print("R rho" , self.rho.R.p)
				# print("D rho" , self.rho.D.p)
				self.RungeKutta(self.rho , i)			
				# print("R rhoU" , self.rhoU.R.p)
				self.RungeKutta(self.rhoU , i)
				# print("R rhoE" , self.rhoE.R.p)
				self.RungeKutta(self.rhoE , i)

			self.rho.residual.append(abs(rho0[3] - self.rho.field[3]))

			# self.rho.print()
			# self.rhoU.print()
			# self.rhoE.print()


			self.correctFields()
			self.correctOutlet()
			self.correctInlet()

			# self.rho.print()
			# self.rhoU.print()
			# self.rhoE.print()

			it = it + 1
			print("-----------------------------")
		
		print("rho residual " , self.rho.residual)

		
	
	
		
	

problem = Euler(4 ,  1e-9 , 1 , 0.1 , 0.3 , 1e5 , 0.9 , 1)

problem.solve()

plt.plot(problem.rho.residual)
plt.yscale('log')
plt.show()

plt.figure()
plt.plot(problem.rho.field)
plt.show()
# problem.rho.printFlux()
# problem.rho.print()
# problem.rhoU.printFlux()
# problem.rhoU.print()
# problem.rhoE.printFlux()
# problem.rhoE.print()
# print()

# problem.computeFluxes()
# problem.rho.printFlux()
# problem.rhoU.printFlux()
# problem.rhoE.printFlux()
# print()

# problem.correctInlet()
# problem.rho.printFlux()
# problem.rho.print()
# problem.rhoU.printFlux()
# problem.rhoU.print()
# problem.rhoE.printFlux()
# problem.rhoE.print()
# print()

# problem.correctOutlet()
# problem.rho.printFlux()
# problem.rho.print()
# problem.rhoU.printFlux()
# problem.rhoU.print()
# problem.rhoE.printFlux()
# problem.rhoE.print()
# print()

# problem.computeFluxes()
# problem.rho.printFlux()
# problem.rhoU.printFlux()
# problem.rhoE.printFlux()


	