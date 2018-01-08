import numpy as np
import matplotlib
import math
import matplotlib.pyplot as plt
from matplotlib import cm
import os

def read_file_and_plot_curve(filename,Nf,Nx,Ny,log_scale,axis,title,name):
	file0 = []
	Nyp   = []
	for n in range(0,Nf):
		file0.append(open(filename+str(n+1)+'.dat', 'r'))
		t0      = []
		counter = 0
		line    = file0[n].readline()
		line    = line.strip()
		array   = line.split()
		t0.append(float(array[0]))
		for line in file0[n]:
			line  = line.strip()
 			array = line.split()
			t0.append(float(array[0]))
			if (t0[counter+1] != t0[counter]):
				Nyp.append(1 + int(counter / Nx))
				file0[n].close
				break
			else:
				counter = counter + 1	
	N  = Nx*Ny
	Np = [] 
	for n in range(0,Nf):
		Np.append(Nx*Nyp[n])
	file    = []
	t       = []
	x       = []
	y       = []
	p       = []
	counter = 0
	if (axis == 1) :
		txt = '_vs_x'
	elif (axis == 2) :
		txt = '_vs_y'
	elif (axis == 3) :
		txt = '_vs_r'
	print('curve '+filename+txt+'.dat :')
	for n in range(0,Nf):
		file.append(open(filename+str(n+1)+'.dat', 'r'))
	for line in file[0]:
		line      = line.strip()
		array     = line.split()
		t.append(float(array[0]))
		y.append(float(array[1]))
		x.append(float(array[2]))
		if isinstance(float(array[3]), float) :
			p.append(float(array[3]))
		else :
			p.append(0.)
		counter = counter + 1
		if (counter % Np[0] == 0):
			print('- at t = '+str(math.floor(100.*t[counter-1])/100))
			M = counter/Np[0]
			for n in range(1,Nf):
				for k in range(0,Np[n]):
					line  = file[n].readline()
					line  = line.strip()
					array = line.split()
					y.append(float(array[1]))
					x.append(float(array[2]))
					if isinstance(float(array[3]), float) :
						p.append(float(array[3]))
					else :
						p.append(math.pow(10.,100.))
			X = np.zeros((Nx,Ny))
			Y = np.zeros((Nx,Ny))
			P = np.zeros((Nx,Ny))
 			XX  = np.zeros(Nx)
 			PPx = np.zeros(Nx)
 			YY  = np.zeros(Ny)
 			PPy = np.zeros(Ny)
			if (Nx <= Ny) :
				Nr = Nx
				NxleNy = 1
			else :
				Nr = Ny
				NxleNy = 0
 			R   = np.zeros(Nr)
 			PPr = np.zeros(Nr)
 			if (log_scale == 1) :
 				Minval = -15.
			for j in range(0,Ny):
				for i in range(0,Nx):
					X[i][j]=x[(M-1)*N+j*Nx+i]
					Y[i][j]=y[(M-1)*N+j*Nx+i]
					if (log_scale == 0) :
						P[i][j]=p[(M-1)*N+j*Nx+i]
					else :
						if (abs(p[(M-1)*N+j*Nx+i]) > math.pow(10.,Minval)) :
							P[i][j]=math.log(np.abs(p[(M-1)*N+j*Nx+i]))/math.log(10.)
						else:
							P[i][j]=Minval
			for j in range(0,Ny):
				for i in range(0,Nx):
					if ((axis == 1) and (j == Ny/2)) :
						XX[i]  = X[i][j]
						PPx[i] = P[i][j]
					if ((axis == 2) and (i == Nx/2)) :
						YY[j]  = Y[i][j]
						PPy[j] = P[i][j]
					if ((axis == 3) and (X[i][j] == Y[i][j]) ) :
						if (NxleNy == 1) :
							r = i
						else :
							r = j
						R0 = np.sqrt((X[Nx/2][1]*X[Nx/2][1])+(Y[1][Ny/2]*Y[1][Ny/2])) 
						R1 = np.sqrt((X[i][j]   *X[i][j]   )+(Y[i][j]   *Y[i][j]   ))
						R[r]   = R1 - R0
						PPr[r] = P[i][j]
			if (axis == 3) :
				aa = np.abs(np.matrix(R).min())
				bb = np.matrix(R).max()
				if ( aa < bb ) :
					Rmax = aa
				else :
					Rmax = bb
			Maxval = np.matrix(P).max()
			Minval = np.matrix(P).min()
			fig=plt.figure()
			font = {'family': 'non-serif',
					'style':  'normal',
        			'color':  'black',
       		 		'weight': 'normal',
       		 		'size': 16,
       				}
			plt.rc('text', usetex=True)
			if (axis == 1) :
				plt.plot(XX,PPx,'k')
				plt.xlim([XX[0],XX[Nx-1]])
				plt.xlabel(r'$x\,(\Delta)$', fontdict=font)
				plt.title(r'$y = $'+str(Y[1][Ny/2])+r'$\Delta,\,t = $'+str(math.floor(100.*t[counter-1])/100.)+r'$/\omega$')
			elif (axis == 2) :
				plt.plot(YY,PPy,'k')
				plt.xlim([YY[0],YY[Ny-1]])
				plt.xlabel(r'$y\,(\Delta)$', fontdict=font)
				plt.title('x = '+str(X[Nx/2][1])+r'$/\Delta,\,t = $'+str(math.floor(100.*t[counter-1])/100.)+r'$/\omega$')
			elif (axis == 3) :
				plt.plot(R,PPr,'k')
				plt.xlim([-Rmax,Rmax])
				plt.xlabel(r'$r\,(\Delta)$', fontdict=font)
				plt.title(r'$\theta = \pi / 4,\,t = $'+str(math.floor(100.*t[counter-1])/100.)+r'$/\omega$')
			if (log_scale == 1) :
				if (axis == 1) :
					Minval = math.floor(np.matrix(PPx).min())-1
					Maxval = math.floor(np.matrix(PPx).max())+1
				elif (axis == 2) :
					Minval = math.floor(np.matrix(PPy).min())-1
					Maxval = math.floor(np.matrix(PPy).max())+1
				elif (axis == 3) :
					Minval = math.floor(np.matrix(PPr).min())-1
					Maxval = math.floor(np.matrix(PPr).max())+1
				if (axis == 3) :
					plt.ylabel(r'$\log_{10}[$'+title+r'$(r,\theta,t) (m_e \Delta {\\omega}^2 / e) ]$', fontdict=font)
				else :
					plt.ylabel(r'$\log_{10}[$'+title+r'$(x,y,t) (m_e \Delta {\\omega}^2 / e) ]$', fontdict=font)
			else :
				if (axis == 3) :
					plt.ylabel(title+'r$(r,\,\theta,\,t)\, (m_e \Delta {\\omega}^2 / e)$', fontdict=font)
				else :
					plt.ylabel(title+'r$(x,\,y,\,t)\, (m_e \Delta {\\omega}^2 / e)$', fontdict=font)
			plt.ylim([Minval,Maxval])
			plt.xticks(fontsize=16)
			plt.yticks(fontsize=16)
			fig.savefig(name+str(M)+'.png')
			plt.close(fig)
	for n in range(0,Nf):
		file[n].close
	return;

def read_file_and_plot_curve_Fourier(filename,Nf,Nx,Ny,log_scale,axis,title,name):
	kxp      = 1 + int( ( (Nx/2) - 1 ) / Nf )
	modesx   = Nx / 4
	modesxpd = min(modesx,kxp)
	Nx       = Nf * modesxpd
	modesy   = Ny / 4
	modesyd  = min(2*modesy-1,Ny)
	Ny       = modesyd
	N        = Nx*Ny
	t        = []
	kx       = []
	ky       = []
	p        = []
	file     = []
	counter  = 0
	print('plot '+filename+'.dat :')
	for n in range(0,Nf):
		file.append(open(filename+str(n+1)+'.dat', 'r'))
	for line in file[0]:
		line      = line.strip()
		array     = line.split()
		kx.append(float(array[1]))
		ky.append(float(array[2]))
		if isinstance(float(array[3]), float) :
			p.append(float(array[3]))
		else :
			p.append(0.)
		t.append(float(array[0]))
		counter = counter + 1
		if (counter % (N/Nf) == 0):
			print('- at t = '+str(math.floor(100.*t[counter-1])/100))
			M = int(counter/(N/Nf))
			for n in range(1,Nf):
				for k in range(0,N/Nf):
					line = file[n].readline()
					line      = line.strip()
					array     = line.split()
					kx.append(float(array[1]))
					ky.append(float(array[2]))
					if isinstance(float(array[3]), float) :
						p.append(float(array[3]))
					else :
						p.append(math.pow(10.,100.))
			if Nf == 1:
				Nx2 = Nx
			else:
				Nx2 = Nx / 2
			Kx = np.zeros((Nx2,Ny-2))
			Ky = np.zeros((Nx2,Ny-2))
			P  = np.zeros((Nx2,Ny-2))
			KkX  = np.zeros(Nx2)
 			PPx = np.zeros(Nx2)
 			KkY  = np.zeros(Ny-2)
 			PPy = np.zeros(Ny-2)
 			if (log_scale == 1) :
 				Minval = -15.
			for j in range(0,Ny-2):
				for i in range(0,Nx2):
					Kx[i][j] = kx[(M-1)*N+i*Ny+j+1]
					Ky[i][j] = ky[(M-1)*N+i*Ny+j+1]
					if (log_scale == 0) :
						P[i][j]=p[ (M-1)*N+i*Ny+j+1]
					else :
						if (abs(p[(M-1)*N+i*Ny+j+1]) > math.pow(10.,Minval)) :
							P[i][j]=math.log(np.abs(p[(M-1)*N+i*Ny+j+1]))/math.log(10.)
						else:
							P[i][j]=Minval
					if ((axis == 1) and (j == 1+(Ny/2))) :
						KkX[i] = Kx[i][j]
						PPx[i] = P[i][j]
					if ((axis == 2) and (i == Nx2/2)) :
						KkY[j] = Ky[i][j]
						PPy[j] = P[i][j]
			Maxval = np.matrix(P).max()
			Minval = np.matrix(P).min()
			fig=plt.figure()
			font = {'family': 'non-serif',
					'style':  'normal',
        			'color':  'black',
       		 		'weight': 'normal',
       		 		'size': 16,
       				}
			plt.rc('text', usetex=True)
			if (axis == 1) :
				plt.plot( KkX,PPx,'k')
				plt.plot(-KkX,PPx,'k')
				plt.xlim([-KkX[Nx2-1],KkX[Nx2-1]])
				plt.xlabel(r'$k_x\,(/\Delta)$', fontdict=font)
				plt.title(r'$k_y = $'+str(math.floor(Ky[1][1+(Ny/2)]))+r'$/\Delta,\,t = $'+str(math.floor(100.*t[counter-1])/100.)+r'$/\omega$')			
			elif (axis == 2) :
				plt.plot(KkY,PPy,'k')
				plt.xlim([KkY[0],KkY[Ny-3]])
				plt.xlabel(r'$k_y\,(/\Delta)$', fontdict=font)
				plt.title(r'$k_x = $'+str(math.floor(Kx[Nx2/2][1]))+r'$/\Delta,\,t = $'+str(math.floor(100.*t[counter-1])/100.)+r'$/\omega$')
			if (log_scale == 1) :
				if (axis == 1) :
					Minval = math.floor(np.matrix(PPx).min())-1
					Maxval = math.floor(np.matrix(PPx).max())+1
				elif (axis == 2) :
					Minval = math.floor(np.matrix(PPy).min())-1
					Maxval = math.floor(np.matrix(PPy).max())+1
				plt.ylabel(r'$\log_{10}[$'+title+r'$(k_x,k_y,t) (m_e \Delta {\\omega}^2 / e) ]$', fontdict=font)
			else :
				plt.ylabel(title+'r$(k_x,\,k_y,\,t)\, (m_e \Delta {\\omega}^2 / e)$', fontdict=font)
			plt.ylim([Minval,Maxval])
			plt.xticks(fontsize=16)
			plt.yticks(fontsize=16)
			fig.savefig(name+str(M)+'.png')
			plt.close(fig)
	for n in range(0,Nf):
		file[n].close
	return;

def TwoD_pcolormesh(filename,Nf,Nx,Ny,log_scale,cmap,title,name):
	file0 = []
	Nyp   = []
	for n in range(0,Nf):
		file0.append(open(filename+str(n+1)+'.dat', 'r'))
		t0      = []
		counter = 0
		line    = file0[n].readline()
		line    = line.strip()
		array   = line.split()
		t0.append(float(array[0]))
		for line in file0[n]:
			line  = line.strip()
 			array = line.split()
			t0.append(float(array[0]))
			if (t0[counter+1] != t0[counter]):
				Nyp.append(1 + int(counter / Nx))
				file0[n].close
				break
			else:
				counter = counter + 1	
	N  = Nx*Ny
	Np = [] 
	for n in range(0,Nf):
		Np.append(Nx*Nyp[n])
	file    = []
	t       = []
	x       = []
	y       = []
	p       = []
	counter = 0
	print('plot '+filename+'.dat :')
	for n in range(0,Nf):
		file.append(open(filename+str(n+1)+'.dat', 'r'))
	for line in file[0]:
		line      = line.strip()
		array     = line.split()
		t.append(float(array[0]))
		y.append(float(array[1]))
		x.append(float(array[2]))
		if isinstance(float(array[3]), float) :
			p.append(float(array[3]))
		else :
			p.append(math.pow(10.,100.))
		counter = counter + 1
		if (counter % Np[0] == 0):
			print('- at t = '+str(math.floor(100.*t[counter-1])/100))
			M = counter/Np[0]
			for n in range(1,Nf):
				for k in range(0,Np[n]):
					line  = file[n].readline()
					line  = line.strip()
					array = line.split()
					y.append(float(array[1]))
					x.append(float(array[2]))
					if isinstance(float(array[3]), float) :
						p.append(float(array[3]))
					else :
						p.append(math.pow(10.,100.))
			X = np.zeros((Nx,Ny))
			Y = np.zeros((Nx,Ny))
			P = np.zeros((Nx,Ny))
 			if log_scale : 
 				Minval = -15.
			for j in range(0,Ny):
				for i in range(0,Nx):
					X[i][j]=x[(M-1)*N+j*Nx+i]
					Y[i][j]=y[(M-1)*N+j*Nx+i]
					if log_scale :
						if (abs(p[(M-1)*N+j*Nx+i]) > math.pow(10.,Minval)) :
							P[i][j]=math.log(np.abs(p[(M-1)*N+j*Nx+i]))/math.log(10.)
						else:
							P[i][j]=Minval
					else :
						P[i][j]=p[(M-1)*N+j*Nx+i]
			cmap = plt.get_cmap(cmap)
			Maxval = np.matrix(P).max()
			if log_scale :
				Maxval = math.floor(Maxval) + 1
			else :
				Minval = np.matrix(P).min()
				if (Maxval == Minval):
					Minval = Minval - 1
					Maxval = Maxval + 1
				elif (filename != 'diag/rho'):
					Maxval = max(abs(Maxval),abs(Minval))
					Minval = - Maxval
			norm = cm.colors.Normalize(vmax=Maxval, vmin=Minval)
			fig=plt.figure()
			plt.rc('text', usetex=True)
			if (log_scale) :
				plt.pcolormesh(X,Y,P,cmap=cmap,vmin=Minval,vmax=Maxval)
			else :
				plt.pcolormesh(X,Y,P,cmap=cmap,norm=norm,vmax=Maxval,vmin=Minval)
# 				plt.pcolormesh(X,Y,P,cmap=cmap,norm=norm,vmax=10,vmin=-10)
# 			X1 = 30.*np.ones(Ny)
# 			Y1 = np.zeros(Ny)
# 			X2 = 98.*np.ones(Ny)
# 			Y2 = np.zeros(Ny)
# 			X3 = np.zeros(Nx)
# 			Y3 = 30.*np.ones(Nx)
# 			X4 = np.zeros(Nx)
# 			Y4 = 98.*np.ones(Nx)
# 			for j in range(0,Ny):
#  				Y1[j]=Y[1][j]
#  			Y2 = Y1
#  			for i in range(0,Nx):
#  				X3[j]=X[i][1]
#  			X4 = X3
# 			plt.plot(X1,Y1,'red',linewidth=2,linestyle='--')
# 			plt.plot(X2,Y2,'red',linewidth=2,linestyle='--')
# 			plt.plot(X3,Y3,'red',linewidth=2,linestyle='--')
# 			plt.plot(X4,Y4,'red',linewidth=2,linestyle='--')
			plt.draw()  			
			cbar=plt.colorbar()
			cbar.ax.tick_params(labelsize=16)
			plt.gca().set_aspect('equal')
			font = {'family': 'non-serif',
					'style':  'normal',
        			'color':  'black',
       		 		'weight': 'normal',
       		 		'size': 16,
       				}
			if log_scale :
				plt.title(r'$\log_{10}( $'+title+r'$)$'+' at t = '+str(math.floor(100.*t[counter-1])/100.)+r'$/\omega$', fontdict=font)
			else :
				plt.title(title+' at t = '+str(math.floor(100.*t[counter-1])/100.)+r'$/\omega$', fontdict=font)
			plt.xticks(fontsize=16)
			plt.xlabel(r'$x\,(\Delta)$', fontdict=font)
			plt.xlim([np.matrix(X).min(),np.matrix(X).max()])
			plt.ylabel(r'$y\,(\Delta)$', fontdict=font)
			plt.yticks(fontsize=16)
			plt.ylim([np.matrix(Y).min(),np.matrix(Y).max()])
			plt.subplots_adjust(bottom=0.1, left = 0.1, right=0.9, top=0.9)
			plt.tight_layout()
			fig.savefig(name+str(M)+'.png',bbox_inches='tight')
			plt.close(fig)
	for n in range(0,Nf):
		file[n].close
	return;

def TwoD_pcolormesh_Fourier_modulus(filename,Nf,Nx,Ny,cmap,title,name):
	kxp      = 1 + int( ( (Nx/2) - 1 ) / Nf )
	modesx   = Nx / 4
	modesxpd = min(modesx,kxp)
	Nx       = Nf * modesxpd
	modesy   = Ny / 4
	modesyd  = min(2*modesy-1,Ny)
	Ny       = modesyd
	N        = Nx*Ny
	kx       = []
	ky       = []
	p        = []
	file     = []
	counter  = 0
	print('plot '+filename+'.dat :')
	for n in range(0,Nf):
		file.append(open(filename+str(n+1)+'.dat', 'r'))
	for line in file[0]:
		line      = line.strip()
		array     = line.split()
		kx.append(float(array[1]))
		ky.append(float(array[2]))
		p.append(float(array[3]))
		time = float(array[0])
		counter = counter + 1
		if (counter % (N/Nf) == 0):
			print('- at t = '+str(math.floor(100.*time)/100))
			M = int(counter/(N/Nf))
			for n in range(1,Nf):
				for k in range(0,N/Nf):
					line = file[n].readline()
					line      = line.strip()
					array     = line.split()
					kx.append(float(array[1]))
					ky.append(float(array[2]))
					p.append(float(array[3]))
			if Nf == 1:
				Nx2 = Nx
			else:
				Nx2 = Nx / 2
			Kx = np.zeros((Nx2,Ny-2))
			Ky = np.zeros((Nx2,Ny-2))
			P  = np.zeros((Nx2,Ny-2))
			Minval = -15
			for j in range(0,Ny-2):
				for i in range(0,Nx2):
					Kx[i][j] = kx[(M-1)*N+i*Ny+j+1]
					Ky[i][j] = ky[(M-1)*N+i*Ny+j+1]
# 					P[i][j]  = p[ (M-1)*N+i*Ny+j+1]
					if (abs(p[(M-1)*N+i*Ny+j+1]) > math.pow(10.,Minval)) :
						P[i][j]=math.log(np.abs(p[(M-1)*N+i*Ny+j+1]))/math.log(10.)
					else:
 						P[i][j]=Minval	
			cmap = plt.get_cmap(cmap)
			Maxval = np.matrix(P).max()
# 			Minval = np.matrix(P).min()
			norm = cm.colors.Normalize(vmax=Maxval, vmin=Minval)
			fig=plt.figure()
			plt.rc('text', usetex=True)
			plt.pcolormesh( Kx,Ky,P,cmap=cmap,vmax=Maxval,vmin=Minval)
			plt.pcolormesh(-Kx,Ky,P,cmap=cmap,vmax=Maxval,vmin=Minval)
# 			plt.pcolormesh( Kx,Ky,P,cmap=cmap,norm=norm,vmax=Maxval,vmin=Minval)
# 			plt.pcolormesh(-Kx,Ky,P,cmap=cmap,norm=norm,vmax=Maxval,vmin=Minval)
			cbar=plt.colorbar()
			cbar.ax.tick_params(labelsize=16)
			plt.gca().set_aspect('equal')
			font = {'family': 'non-serif',
					'style':  'normal',
        			'color':  'black',
       			 	'weight': 'normal',
       		 		'size': 16,
       				}
			plt.title(title+' at t = '+str(math.floor(100.*time)/100.)+r'$/\omega$', fontdict=font)
			plt.xticks(fontsize=16)
			plt.xlabel(r'$k_x\,(/\Delta)$', fontdict=font)
			plt.xlim([-np.matrix(Kx).max(),np.matrix(Kx).max()])
			plt.ylabel(r'$k_y\,(/\Delta)$', fontdict=font)
			plt.yticks(fontsize=16)
			plt.ylim([np.matrix(Ky).min(),np.matrix(Ky).max()])
			fig.savefig(name+str(M)+'.png',bbox_inches='tight')
			plt.close(fig)
	for n in range(0,Nf):
		file[n].close
	return;
	
def TwoD_pcolormesh_Fourier_argument(filename,Nf,Nx,Ny,cmap,title,name):
	kxp      = 1 + int( ( (Nx/2) - 1 ) / Nf )
	modesx   = Nx / 4
	modesxpd = min(modesx,kxp)
	Nx       = Nf * modesxpd
	modesy   = Ny / 4
	modesyd  = min(2*modesy-1,Ny)
	Ny       = modesyd
	N        = Nx*Ny
	kx       = []
	ky       = []
	p        = []
	file     = []
	counter  = 0
	print('plot '+filename+'.dat :')
	for n in range(0,Nf):
		file.append(open(filename+str(n+1)+'.dat', 'r'))
	for line in file[0]:
		line      = line.strip()
		array     = line.split()
		kx.append(float(array[1]))
		ky.append(float(array[2]))
		p.append(float(array[3]))
		time = float(array[0])
		counter = counter + 1
		if (counter % (N/Nf) == 0):
			print('- at t = '+str(math.floor(100.*time)/100))
			M = int(counter/(N/Nf))
			for n in range(1,Nf):
				for k in range(0,N/Nf):
					line = file[n].readline()
					line      = line.strip()
					array     = line.split()
					kx.append(float(array[1]))
					ky.append(float(array[2]))
					p.append(float(array[3]))
			if Nf == 1:
				Nx2 = Nx
			else:
				Nx2 = Nx / 2
			Kx = np.zeros((Nx2,Ny-2))
			Ky = np.zeros((Nx2,Ny-2))
			P = np.zeros((Nx2,Ny-2))
			for j in range(0,Ny-2):
				for i in range(0,Nx2):
					Kx[i][j] = kx[(M-1)*N+i*Ny+j+1]
					Ky[i][j] = ky[(M-1)*N+i*Ny+j+1]
					P[i][j]  = p[(M-1)*N+i*Ny+j+1]*(360./6.28318530717959)
			cmap = plt.get_cmap(cmap)
			Maxval = np.matrix(P).max()
			Minval = np.matrix(P).min()
			norm = cm.colors.Normalize(vmax=Maxval, vmin=Minval)
			fig=plt.figure()
			plt.rc('text', usetex=True)
			plt.pcolormesh( Kx,Ky,P,cmap=cmap,norm=norm,vmax=Maxval,vmin=Minval)
			plt.pcolormesh(-Kx,Ky,P,cmap=cmap,norm=norm,vmax=Maxval,vmin=Minval)
			cbar=plt.colorbar()
			cbar.ax.tick_params(labelsize=16)
			plt.gca().set_aspect('equal')
			font = {'family': 'non-serif',
					'style':  'normal',
        			'color':  'black',
       			 	'weight': 'normal',
       		 		'size': 16,
       				}
			plt.title(title+' at t = '+str(math.floor(100.*time)/100.)+r'$/\omega$', fontdict=font)
			plt.xticks(fontsize=16)
			plt.xlabel(r'$k_x\,(/\Delta)$', fontdict=font)
			plt.xlim([-np.matrix(Kx).max(),np.matrix(Kx).max()])
			plt.ylabel(r'$k_y\,(/\Delta)$', fontdict=font)
			plt.yticks(fontsize=16)
			plt.ylim([np.matrix(Ky).min(),np.matrix(Ky).max()])
			fig.savefig(name+str(M)+'.png',bbox_inches='tight')
			plt.close(fig)
	for n in range(0,Nf):
		file[n].close
	return;

def part_plot(dir,filename,Nf,name):
	print('plot '+filename+'.dat :')
	print('1) Find the number of time iterations to be plotted')
	print('...')
	file0 = open(filename+str(1)+'.dat', 'r')
	t0      = []
	counter0 = 0
	line0    = file0.readline()
	line0    = line0.strip()
	array0   = line0.split()
	t0.append(float(array0[0]))
	Nit = 1
	for line0 in file0:
		line0  = line0.strip()
 		array0 = line0.split()
		t0.append(float(array0[0]))
		if (t0[counter0+1] != t0[counter0]):
			Nit += 1
		counter0 +=1
	file0.close()
	file1 = []
	Npp   = np.zeros((Nit,Nf))
	print('2) Find the number of particles per MPI node to plot at each time step')
	print('...')
	for n in range(0,Nf):
		file1.append(open(filename+str(n+1)+'.dat', 'r'))
	t1       = []
	counter1 = 0
	line1    = file1[0].readline()
	line1    = line1.strip()
	array1   = line1.split()
	t1.append(float(array1[0]))
	m = 0
	counter2 = []
	time = np.zeros(Nit)
	for line1 in file1[0]:
		line1  = line1.strip()
 		array1 = line1.split()
		t1.append(float(array1[0]))
		if (t1[counter1+1] != t1[counter1]):
			Npp[m][0] = (1 + counter1)/(m+1)
			for n in range(1,Nf):
				t2     = []
				line2  = file1[n].readline()
				line2  = line2.strip()
				array2 = line2.split()
				counter2.append(1)
				counter2[n-1] = 1
				t2.append(float(array2[0]))
				line2  = file1[n].readline()
				line2  = line2.strip()
				array2 = line2.split()
				t2.append(float(array2[0]))
				while (t2[counter2[n-1]] == t2[counter2[n-1]-1]):
					line2  = file1[n].readline()
					line2  = line2.strip()
					array2 = line2.split()
					t2.append(float(array2[0]))
					counter2[n-1] += 1
				if (m == 0):
					Npp[m][n] = counter2[n-1]
				else:
					Npp[m][n] = 1 + counter2[n-1]
			time[m] = t1[counter1]
			m += 1
		counter1 +=1
	Npp[m][0] = (1 + counter1)/(m+1)
	time[m] = t1[counter1]
	for n in range(1,Nf):
		Npp[m][n] = 1 + counter2[n-1]
	for n in range(0,Nf):
		file1[n].close()
	#print(Npp)
	print('3) Plot')
	file3 = []
	for n in range(0,Nf):
		file3.append(open(filename+str(n+1)+'.dat', 'r'))
	for m in range(0,Nit):
		print('- at t = '+str(math.floor(100.*time[m])/100))
		Npp_m = np.zeros(Nf)
		Npp_m = Npp[m][0:Nf]
		Npp_tot = np.sum(Npp_m[0:Nf])
		Npp_tot = int(Npp_tot)
		x = np.zeros(Npp_tot)
		y = np.zeros(Npp_tot)
		p = np.zeros(Npp_tot)
		K = 0
		for n in range(0,Nf):
			for k in range(0,int(Npp[m][n])):
				line3  = file3[n].readline()
				line3  = line3.strip()
 				array3 = line3.split()
				x[K+k]=float(array3[1])
				y[K+k]=float(array3[2])
				p[K+k]=float(array3[3])
			K = K + int(Npp[m][n])
		fig = plt.figure()
		plt.rc('text', usetex=True)
		font = {'family': 'non-serif',
				'style':  'normal',
        		'color':  'black',
       			 'weight': 'normal',
       			 'size': 16,
       			}
		ax = fig.add_subplot(111, projection='3d')
		cm = plt.get_cmap("jet")
		col = [cm(float(i)/len(p)) for i in xrange(len(p))]
		ax.plot(x, y, p,'bo', markersize=2.,alpha=0.25)
		#ax.scatter(x, y, p, s=20, c=col, marker='o', edgecolors=None,alpha=0.125)
		ax.set_xlim3d([np.amin(x),np.amax(x)])
		ax.set_xlabel(r'$x\,(\Delta)$', fontdict=font)
		ax.set_ylim3d([np.amin(y),np.amax(y)])
		ax.set_ylabel(r'$y\,(\Delta)$', fontdict=font)
		ax.set_zlim3d([np.amin(p),np.amax(p)])
		#ax.set_zlim3d([-12,12])
		if (dir == 1):
			ax.set_zlabel(r'$v_x \,(\Delta \\omega)$', fontdict=font)
			#ax.view_init(0,-90)
		elif (dir == 2):
			ax.set_zlabel(r'$v_y \,(\Delta \\omega$)', fontdict=font)
			#ax.view_init(0,0)
		else :
			ax.set_zlabel(r'$v_z \,(\Delta \\omega$)', fontdict=font)
		plt.title('phase-space at '+str(math.floor(100.*time[m])/100)+r'$\,(/\omega)$', fontdict=font)
		fig.savefig(name+str(m+1),bbox_inches='tight')
 		plt.close(fig)
	for n in range(0,Nf):
		file3[n].close
	return;

def read_file_and_define_first_col(filename,time):
	file = open(filename, 'r')
	for line in file:
		line      = line.strip()
		array     = line.split()
		time.append(float(array[0]))
	file.close()   
	return;

def read_file_and_define_second_col(filename,vector):
	file = open(filename, 'r')
	for line in file:
		line      = line.strip()
		array     = line.split()
		F = float(array[1])
		if (abs(F) > 1.e-15) :
			vector.append(F)
		else:
 			vector.append(1.e-15)
	file.close()   
	return;