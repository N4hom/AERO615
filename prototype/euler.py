import numpy as np

class variable(object):
	"""docstring for field"""
	def __init__(self, name : str  , Nx, Ny , top, bottom, left, right ):
		
		self.name = name

		self.Nx = Nx
		self.Ny = Ny

		
		# Define coefficient matrix for the variables to be solved
		self.a = coeff(np.zeros((Nx,Ny)) , np.zeros((Nx,Ny)) , np.zeros((Nx,Ny)) , np.zeros((Nx,Ny)), np.zeros((Nx,Ny)), np.zeros((Nx,Ny)))
		
		self.iMax = Nx + 1
		self.jMax = Ny + 1
		self.field = np.zeros((Nx+2,Ny+2))


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


		
		self.iMax = Nx + 1
		# self.jMax = Ny + 1
		self.field = np.ones(Nx+2)*value
		self.flux_f = np.zeros(Nx+2)
		self.flux_g = np.zeros(Nx+2)


		
		self.left = left
		self.right = right
		
	def __getitem__(self, i):
		return self.field[i]

	def __setitem__(self, i, val):
		self.field[i] = val

	def print(self):
		print(self.name ,self.field)
		return 

	def interpolateLeft(self , i ):
		return 0.5 * (self.field[i] + self.field[i - 1])

	def interpolateRight(self , i ):
		return 0.5 * (self.field[i] + self.field[i + 1])

	def computeResidualij(self , i):
		self.R.w[i] = - self.interpolateLeft(i) * 0.1  # face length is hardcoded
		self.R.e[i] =   self.interpolateRight(i) * 0.1

		self.R.p[i] = self.R.w[i] + self.R.e[i]

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

		self.pInf = pInf
		self.rhoInf = rhoInf
		self.pRatio = pRatio
		self.M = Mach
		self.cInf = np.sqrt(1.4 * self.pInf / self.rhoInf)
		self.uInf = self.M * self.cInf
		self.gamma = gamma
		self.gamma_1 = gamma - 1


		self.rho = variable1D  ("rho"   , Nx , 0.0 , 0.0 , self.rhoInf * 0.9)
		self.rhoU = variable1D ("rhoU"  , Nx , 0.  , 0.  , self.rhoInf * self.uInf )  
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

		Riem1[0] = uInf + 2/gamma_1 * cInf
		Riem2[0] = u[1] - 2/gamma_1 * c[1]  # Riem2[0] = Riem2[1]  

		u[0] = 0.5 * (Riem1[0] + Riem2[0])
		c[0] = 0.25 * gamma_1 * (Riem1[0] - Riem2[0])
		M[0] = u[0]/c[0]
		p[0] = pInf/(1 + 0.5 * gamma_1 * M[0]**2)**(gamma/gamma_1)
		rho[0] = gamma * p[0]/c[0]**2
		rhoE[0] = p[0]/gamma_1 + 0.5 * rho[0] * u[0]**2
		rhoU[0] = rho[0] * u[0]


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

		print("Correct outlet")
		p[-2] = pInf*self.pRatio
		p[-1] = p[-2]

		rho[-2] = 2*rho[-3] - rho[-1]
		rho[-1] = rho[-2]

		rhoU[-2] = 2*rhoU[-3] - rhoU[-1]
		rhoU[-1] = rhoU[-2]

		rhoE[-2] = p[-1]/gamma_1 + 0.5 * rho[-2] * u[-2]**2
		rhoE[-1] = rhoE[-2] 

		c[-2] = np.sqrt(gamma * p[-2] / rho[-2])
		c[-1] = np.sqrt(gamma * p[-1] / rho[-1])

		u[-2] = rhoU[-2]/rho[-2]
		u[-1] = u[-2]


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
		u.field[1:-2] = rhoU.field[1:-2]/rho.field[1:-2]
		p.field[1:-2] = gamma_1 * (rhoE.field[1:-2] - 0.5 * rho.field[1:-2] * u.field[1:-2] ** 2)
		c.field[1:-2] = np.sqrt(gamma * p.field[1:-2] / rho.field[1:-2])

		self.p.print()
		self.c.print()

	def computeFluxes(self):
		
		print("Compute fluxes")
		self.rho.computeFlux(self.u.field , None)
		self.rhoU.computeFlux(self.u.field , self.p.field)
		self.rhoE.computeFlux(self.u.field , self.p.field)
		self.u.print()
		print(self.rho.flux_f)
		print(self.rhoU.flux_f)
		print(self.rhoE.flux_f)
	
	
	
		
	

problem = Euler(4 ,  1e-9 , 1 , 0.1 , 0.3 , 1e5 , 0.9 , 1)

print(problem.rho.R.w)

problem.rho.print()
problem.rhoU.print()
problem.rhoE.print()
problem.Riem1.print()
print()
print()

problem.correctInlet()	
problem.rho.print()
problem.rhoU.print()
problem.rhoE.print()
problem.Riem1.print()
print()

problem.correctOutlet()	
problem.rho.print()
problem.rhoU.print()
problem.rhoE.print()
problem.Riem1.print()
print()

problem.correctFields()
problem.computeFluxes()
print()



	