import os
import numpy as np
import read_and_plot
import matplotlib.pyplot as plt
import sys

folder_name = 'plots'
try:
    N_nodes=int(raw_input('How many MPI nodes have been used? '))
except ValueError:
    print "Not a number"
    sys.exit()

try:
    Delta_x=float(raw_input('Delta_x = '))
except ValueError:
    print "Not a number"
    sys.exit()

try:
    Delta_y=float(raw_input('Delta_y = '))
except ValueError:
    print "Not a number"
    sys.exit()

try:
    indx=int(raw_input('indx = '))
except ValueError:
    print "Not a number"
    sys.exit()
    
try:
    indy=int(raw_input('indy = '))
except ValueError:
    print "Not a number"
    sys.exit()

Nx = 2**indx
Ny = 2**indy

try:
    test=int(raw_input('Mobile ions (1 if true or 0 if false)? '))
except ValueError:
    print "Not a number"
    sys.exit()
    
try:
    test2=int(raw_input('Do you want to plot Fourier transforms? (1 : yes or 0 : no)? '))
except ValueError:
    print "Not a number"
    sys.exit()

if (test == 0):
	movion = False
else :
	movion = True

if (test2 == 0):
	spectral = False
else :
	spectral = True

dir= os.path.dirname(folder_name+'/')
if not os.path.exists(dir):
    os.mkdir(dir)

#######################
# What I want to plot #
#######################

# plot_histograms = True
plot_histograms = False

# plot_energies = False
plot_energies = True

# plot_densities = False
plot_densities = True

# densities_log_scale = True
densities_log_scale = False

# plot_fields = False
plot_fields = True

# fields_log_scale = True
fields_log_scale = False

plot_Poynting = False
# plot_Poynting = True

# Poynting_log_scale = True
Poynting_log_scale = False

pml_tests_plot = False
# pml_tests_plot = True

phase_spaces_plot = False
# phase_spaces_plot = True

########################

results_folder = 'results'

if plot_histograms :

####################
#  Histogram test  #
####################

	print('Histogram test :')
	print('-------------')

	dir= os.path.dirname(folder_name+'/histograms/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	dir= os.path.dirname(folder_name+'/histograms/ele/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	if ( movion ):
		dir= os.path.dirname(folder_name+'/histograms/ion/')
		if not os.path.exists(dir):
			os.mkdir(dir)

	for k in range(0,N_nodes) :
		v = []
		F = []
		read_and_plot.read_file_and_define_first_col(results_folder+'/histogram_ele'+str(k+1)+'.dat',v)
		read_and_plot.read_file_and_define_second_col(results_folder+'/histogram_ele'+str(k+1)+'.dat',F)
		n = len(v)
		dv = v[2]-v[1]
		Checksum = 0.
		for i in range(0,n) :
			Checksum = Checksum + (F[i]*dv)
		Ff = np.zeros(n)
		for i in range(0,n) :
			Ff[i] = F[i] / Checksum
		fig=plt.figure()
		plt.rc('text', usetex=True)
		plt.plot(v, Ff,'red',linewidth=2,label=r'$f_e(p) = \frac{\mu_e}{K_2(\mu_e)} \frac{p^2}{(m_ec)^3} \exp{(-\mu_e \gamma(p))},\,\mu_e=\frac{m_e c^2}{k_B T_e}$')
		font = {'family': 'non-serif',
				'style':  'normal',
				'color':  'black',
				'weight': 'normal',
				'size': 16,
				}
		leg = plt.legend(loc='upper right',fontsize=16, fancybox=True)
		leg.get_frame().set_alpha(0.5)
		plt.xticks(fontsize=16)
		plt.xlabel(r'$p\,(m_e \Delta \omega)$', fontdict=font)
		plt.xlim([v[0],v[n-1]])
		plt.ylabel(r'$f_e \,(1 / m_e \Delta \omega)$', fontdict=font)
		plt.yticks(fontsize=16)
		plt.ylim([0,1.1*np.amax(Ff)])
		plt.title('Electron momentum amplitude distribution in the spatial cell (1,'+str(int((k*Ny/N_nodes)+3))+')')
		fig.savefig(folder_name+'/histograms/ele/histogram_ele'+str(int(k+1))+'.png',bbox_inches='tight')
		plt.close(fig)

		if ( movion ):
			vi = []
			Fi = []
			read_and_plot.read_file_and_define_first_col(results_folder+'/histogram_ion'+str(k+1)+'.dat',vi)
			read_and_plot.read_file_and_define_second_col(results_folder+'/histogram_ion'+str(k+1)+'.dat',Fi)	
			n = len(vi)
			dvi = vi[2]-vi[1]
			Checksumi = 0.
			for i in range(0,n) :
				Checksumi = Checksumi + (Fi[i]*dvi)
			Ffi = np.zeros(n)
			for i in range(0,n):
				Ffi[i] = Fi[i] / Checksumi
			fig=plt.figure()
			plt.rc('text', usetex=True)
			plt.plot(vi, Ffi,'green',linewidth=2,label=r'$f_i(p) = \frac{\mu_i}{K_2(\mu_i)} \frac{p^2}{(m_ic)^3} \exp{(-\mu_i \gamma(p))},\,\mu_i=\frac{m_i c^2}{k_B T_i}$')
			leg = plt.legend(loc='upper right',fontsize=16, fancybox=True)
			leg.get_frame().set_alpha(0.5)
			font = {'family': 'non-serif',
					'style':  'normal',
					'color':  'black',
					'weight': 'normal',
					'size': 16
					}
			plt.xticks(fontsize=16)
			plt.xlabel(r'$p\,(m_i \Delta \omega)$', fontdict=font)
			plt.xlim([vi[0],vi[n-1]])
			plt.ylabel(r'$f_i \,(1 / m_i \Delta \omega)$', fontdict=font)
			plt.yticks(fontsize=16)
			plt.ylim([0,1.1*np.amax(Ffi)])
			plt.title('Ion momentum amplitude distribution in the spatial cell (3,'+str(int(k*(Ny/N_nodes)+3))+')')
			fig.savefig(str(folder_name+'/histograms/ion/histogram_ion'+str(k+1)+'.png'),bbox_inches='tight')
			plt.close(fig)

##############
#  Energies  #
##############
   
if plot_energies :

	print('Energy plot :')
	print('-------------')

	t = []
	read_and_plot.read_file_and_define_first_col(results_folder+'/Ues.dat',t)
	Ues = []
	read_and_plot.read_file_and_define_second_col(results_folder+'/Ues.dat',Ues)
	Uk = []
	read_and_plot.read_file_and_define_second_col(results_folder+'/Uk.dat',Uk)
	if movion :
		Uke = []
		read_and_plot.read_file_and_define_second_col(results_folder+'/Uke.dat',Uke)
		Uki = []
		read_and_plot.read_file_and_define_second_col(results_folder+'/Uki.dat',Uki)
	Usrc = []
 	read_and_plot.read_file_and_define_second_col(results_folder+'/Usrc.dat',Usrc)
	Uel = []
	read_and_plot.read_file_and_define_second_col(results_folder+'/Uel.dat',Uel)
	Uma = []
	read_and_plot.read_file_and_define_second_col(results_folder+'/Uma.dat',Uma)
	Uesc = []
	read_and_plot.read_file_and_define_second_col(results_folder+'/Uesc.dat',Uesc)
	Uel_pml = []
	read_and_plot.read_file_and_define_second_col(results_folder+'/Uel_pml.dat',Uel_pml)
	Uma_pml = []
	read_and_plot.read_file_and_define_second_col(results_folder+'/Uma_pml.dat',Uma_pml)
	Uel_dumped = []
	read_and_plot.read_file_and_define_second_col(results_folder+'/Uel_dumped.dat',Uel_dumped)
	Uma_dumped = []
	read_and_plot.read_file_and_define_second_col(results_folder+'/Uma_dumped.dat',Uma_dumped)
	n           = len(t)
	Utot        = np.zeros(n)
	Utot_pml    = np.zeros(n)
	Utot_dumped = np.zeros(n)
	for i in range(0,n):
		Utot[i]        = Uel[i] + Uma[i] + Uk[i]
		if len(Uel_pml) == n :
			Utot_pml[i]    = Uel_pml[i] + Uma_pml[i]
		if len(Uel_dumped) == n :
			Utot_dumped[i] = Uel_dumped[i] + Uma_dumped[i]
	if (Uesc[n-1] != 0.) :
		error_pml = (Utot_pml[n-1] + Uesc[n-1] - Utot_dumped[n-1] ) / Uesc[n-1]
	if (Usrc[n-1] != 0.) :
		error_sim = (Utot[n-1] - Uesc[n-1] - Usrc[n-1] ) / Usrc[n-1]
		error     = (Utot[n-1] + Utot_pml[n-1] - (Usrc[n-1] + Utot_dumped[n-1]) ) / Usrc[n-1]
	fig=plt.figure()
	plt.rc('text', usetex=True)
	plt.semilogy(t, Ues,'cyan',linewidth=2,label=r'$U_\mathrm{es}$')
# 	plt.semilogy(t, Usrc,'red',linewidth=2,label='time integrated '+r'$U_\mathrm{src}$')
	plt.semilogy(t, Uel,'blue',linewidth=2,label='instantaneous '+r'$U_\mathrm{E}$')
	plt.semilogy(t, Uma,'magenta',linewidth=2,label='instantaneous '+r'$U_\mathrm{B}$')
	plt.semilogy(t, Utot,'black',linewidth=2,label=r'instantaneous $U_\mathrm{E,SIM} + U_\mathrm{B,SIM}$')
# 	plt.semilogy(t, Uesc,'magenta',linewidth=2,label='time integrated '+r'$U_\mathrm{esc}$')
# 	plt.semilogy(t, Uel_pml,'blue',linestyle='--',linewidth=2,label='instantaneous '+r'$U_\mathrm{E,PML}$')
# 	plt.semilogy(t, Uma_pml,'magenta',linestyle='--',linewidth=2,label='instantaneous '+r'$U_\mathrm{B,PML}$')
# 	if len(Uel_pml) == n :
# 		plt.semilogy(t, Utot_pml,'blue',linewidth=2,label='instantaneous '+r'$U_\mathrm{E,PML} + U_\mathrm{B,PML}$')
# 	plt.semilogy(t, Uel_dumped,'cyan',linewidth=2,label='time integrated '+r'$U_\mathrm{dumped,E}$')
# 	plt.semilogy(t, Uma_dumped,'orange',linewidth=2,label='time integrated '+r'$U_\mathrm{dumped,B}$')
# 	if len(Uel_dumped) == n :
# 		plt.semilogy(t, Utot_dumped,'green',linewidth=2,label='time integrated '+r'$U_\mathrm{dumped}$')
	plt.semilogy(t, Uk,'green',linewidth=2,label='time integrated '+r'$U_\mathrm{kin}$')
	if movion :
		plt.semilogy(t, Uke,'green',linewidth=2,label='time integrated '+r'$U_\mathrm{kin,e}$')
		plt.semilogy(t, Uki,'orange',linewidth=2,label='time integrated '+r'$U_\mathrm{kin,i}$')
	leg = plt.legend(loc='lower right',fontsize=16, fancybox=True)
	leg.get_frame().set_alpha(0.5)
	font = {'family': 'non-serif',
			'style':  'normal',
			'color':  'black',
			'weight': 'normal',
			'size': 16,
			}
# 	if len(Uel_dumped) == n :
# 		plt.title(r'$\varepsilon_\mathrm{SIM} = $'+str(np.floor(100*error_sim)/100)+r'$\,\%,\,$'+r'$\varepsilon_\mathrm{PML} = $'+str(np.floor(100*error_pml)/100)+r'$\,\%,\,$'+r'$\varepsilon = $'+str(np.floor(100*error)/100)+r'$\,\%$', fontdict=font)
# 	else :
#  		plt.title(r'$\varepsilon = $'+str(np.floor(100*error)/100)+r'$\,\%$', fontdict=font)
	plt.xticks(fontsize=16)
	plt.xlabel(r'$t\,(/\omega)$', fontdict=font)
	plt.xlim([t[0],t[n-1]])
	plt.ylabel(r'$\mathrm{Energy}\,(m_e \Delta^2 \omega^2)$', fontdict=font)
	plt.yticks(fontsize=16)
	plt.ylim([0,1.1*np.amax(Utot)])
	fig.savefig(folder_name+'/energy.png',bbox_inches='tight')
	plt.close(fig)

#############
# Densities #
#############

if plot_densities :

	print('Density plots :')
	print('---------------')

	dir= os.path.dirname(folder_name+"/rho/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	dir= os.path.dirname(folder_name+"/rho/rho/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.TwoD_pcolormesh(results_folder+'/rho',
							N_nodes,Nx,Ny,densities_log_scale,
							'nipy_spectral',
							r'$\mathbf{\rho}\,(e/\Delta^3)$',
							folder_name+'/rho/rho/rho')

	if movion :
		dir= os.path.dirname(folder_name+"/rho/rho_ele/")
		if not os.path.exists(dir):
			os.mkdir(dir)
	
		read_and_plot.TwoD_pcolormesh(results_folder+'/rho_ele',
							N_nodes,Nx,Ny,densities_log_scale,
							'nipy_spectral',
							r'$\mathbf{\rho}_e\,(e/\Delta^3)$',
							folder_name+'/rho/rho_ele/rho_ele')
		dir= os.path.dirname(folder_name+"/rho/rho_ion/")
		if not os.path.exists(dir):
			os.mkdir(dir)
	
		read_and_plot.TwoD_pcolormesh(results_folder+'/rho_ion',
							N_nodes,Nx,Ny,densities_log_scale,
							'nipy_spectral',
							r'$\mathbf{\rho}_i\,(e/\Delta^3)$',
							folder_name+'/rho/rho_ion/rho_ion')

	if spectral :
		dir= os.path.dirname(folder_name+"/rho/rho_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/rho/rho_Fourier/modulus/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_modulus(results_folder+'/rho_Fourier',
								N_nodes,Nx,Ny,
								'hot',
								r'$|\widehat{\mathbf{\rho}}|\,(e/\Delta^3)$',
								folder_name+'/rho/rho_Fourier/modulus/rho_F_mod')

		dir= os.path.dirname(folder_name+"/rho/rho_Fourier/argument/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.TwoD_pcolormesh_Fourier_argument(results_folder+'/rho_Fourier_arg',
								N_nodes,Nx,Ny,
								'jet',
								r'$\mathrm{arg}\left (\widehat{\mathbf{\rho}}\right )\,(^o)$',
								folder_name+'/rho/rho_Fourier/argument/rho_F_arg')


	dir= os.path.dirname(folder_name+"/j/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	dir= os.path.dirname(folder_name+"/j/j/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	dir= os.path.dirname(folder_name+"/j/j/jx/")
	if not os.path.exists(dir):
		os.mkdir(dir)


	read_and_plot.TwoD_pcolormesh(results_folder+'/jx',
							N_nodes,Nx,Ny,densities_log_scale,
							'Spectral',
							r'$\mathbf{j_x}\,(e \omega /\Delta^2)$',
							folder_name+'/j/j/jx/jx')

	dir= os.path.dirname(folder_name+"/j/j/jy/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.TwoD_pcolormesh(results_folder+'/jy',
							N_nodes,Nx,Ny,densities_log_scale,
							'Spectral',
							r'$\mathbf{j_y}\,(e \omega /\Delta^2)$',
							folder_name+'/j/j/jy/jy')

	dir= os.path.dirname(folder_name+"/j/j/jz/")
	if not os.path.exists(dir):
		os.mkdir(dir)
						
	read_and_plot.TwoD_pcolormesh(results_folder+'/jz',
							N_nodes,Nx,Ny,densities_log_scale,
							'nipy_spectral',#'Spectral',
							r'$\mathbf{j_z}\,(e \omega /\Delta^2)$',
							folder_name+'/j/j/jz/jz')

	if movion :
		dir= os.path.dirname(folder_name+"/j/je/")
		if not os.path.exists(dir):
			os.mkdir(dir)
	
		dir= os.path.dirname(folder_name+"/j/je/jex/")
		if not os.path.exists(dir):
			os.mkdir(dir)


		read_and_plot.TwoD_pcolormesh(results_folder+'/jex',
								N_nodes,Nx,Ny,densities_log_scale,
								'Spectral',
								r'$\mathbf{j_{e,x}}\,(e \omega /\Delta^2)$',
								folder_name+'/j/je/jex/jex')

		dir= os.path.dirname(folder_name+"/j/je/jey/")
		if not os.path.exists(dir):
			os.mkdir(dir)
	
		read_and_plot.TwoD_pcolormesh(results_folder+'/jey',
								N_nodes,Nx,Ny,densities_log_scale,
								'Spectral',
								r'$\mathbf{j_{e,y}}\,(e \omega /\Delta^2)$',
								folder_name+'/j/je/jey/jey')

		dir= os.path.dirname(folder_name+"/j/je/jez/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh(results_folder+'/jez',
								N_nodes,Nx,Ny,densities_log_scale,
								'Spectral',
								r'$\mathbf{j_{e,z}}\,(e \omega /\Delta^2)$',
								folder_name+'/j/je/jez/jez')
							
		dir= os.path.dirname(folder_name+"/j/ji/")
		if not os.path.exists(dir):
			os.mkdir(dir)
	
		dir= os.path.dirname(folder_name+"/j/ji/jix/")
		if not os.path.exists(dir):
			os.mkdir(dir)


		read_and_plot.TwoD_pcolormesh(results_folder+'/jix',
								N_nodes,Nx,Ny,densities_log_scale,
								'Spectral',
								r'$\mathbf{j_{i,x}}\,(e \omega /\Delta^2)$',
								folder_name+'/j/ji/jix/jix')

		dir= os.path.dirname(folder_name+"/j/ji/jiy/")
		if not os.path.exists(dir):
			os.mkdir(dir)
	
		read_and_plot.TwoD_pcolormesh(results_folder+'/jiy',
								N_nodes,Nx,Ny,densities_log_scale,
								'Spectral',
								r'$\mathbf{j_{i,y}}\,(e \omega /\Delta^2)$',
								folder_name+'/j/ji/jiy/jiy')

		dir= os.path.dirname(folder_name+"/j/ji/jiz/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh(results_folder+'/jiz',
								N_nodes,Nx,Ny,densities_log_scale,
								'Spectral',
								r'$\mathbf{j_{i,z}}\,(e \omega /\Delta^2)$',
								folder_name+'/j/ji/jiz/jiz')

	if spectral :
		dir= os.path.dirname(folder_name+"/j/j_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/j/j_Fourier/jx_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/j/j_Fourier/jx_Fourier/modulus/")
		if not os.path.exists(dir):
			os.mkdir(dir)
					
		read_and_plot.TwoD_pcolormesh_Fourier_modulus(results_folder+'/jx_Fourier',
								N_nodes,Nx,Ny,
								'hot',
								r'$|\widehat{\mathbf{j}}_x|\,(e \omega /\Delta^2)$',
								folder_name+'/j/j_Fourier/jx_Fourier/modulus/jx_F_mod')

		dir= os.path.dirname(folder_name+"/j/j_Fourier/jx_Fourier/argument/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.TwoD_pcolormesh_Fourier_argument(results_folder+'/jx_Fourier_arg',
								N_nodes,Nx,Ny,
								'jet',
								r'$\mathrm{arg}\left (\widehat{\mathbf{j}}_x\right)\,(^o)$',
								folder_name+'/j/j_Fourier/jx_Fourier/argument/jx_F_arg')

		dir= os.path.dirname(folder_name+"/j/j_Fourier/jy_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/j/j_Fourier/jy_Fourier/modulus/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_modulus(results_folder+'/jy_Fourier',
								N_nodes,Nx,Ny,
								'hot',
								r'$|\widehat{\mathbf{j}}_y|\,(e \omega /\Delta^2)$',
								folder_name+'/j/j_Fourier/jy_Fourier/modulus/jy_F_mod')

		dir= os.path.dirname(folder_name+"/j/j_Fourier/jy_Fourier/argument/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.TwoD_pcolormesh_Fourier_argument(results_folder+'/jy_Fourier_arg',
								N_nodes,Nx,Ny,
								'jet',
								r'$\mathrm{arg}\left (\widehat{\mathbf{j}}_y\right)\,(^o)$',
								folder_name+'/j/j_Fourier/jy_Fourier/argument/jy_F_arg')

		dir= os.path.dirname(folder_name+"/j/j_Fourier/jz_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/j/j_Fourier/jz_Fourier/modulus/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_modulus(results_folder+'/jz_Fourier',
								N_nodes,Nx,Ny,
								'hot',
								r'$|\widehat{\mathbf{j}}_z|\,(e \omega /\Delta^2)$',
								folder_name+'/j/j_Fourier/jz_Fourier/modulus/jz_F_mod')

		dir= os.path.dirname(folder_name+"/j/j_Fourier/jz_Fourier/argument/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.TwoD_pcolormesh_Fourier_argument(results_folder+'/jz_Fourier_arg',
								N_nodes,Nx,Ny,
								'jet',
								r'$\mathrm{arg}\left (\widehat{\mathbf{j}}_z\right)\,(^o)$',
								folder_name+'/j/j_Fourier/jz_Fourier/argument/jz_F_arg')						

##########################
# Electromagnetic fields #
##########################

if plot_fields :

	print('Electromagnetic fields plots')
	print('----------------------------')

	if fields_log_scale :
		fields_plots_color = 'nipy_spectral'
	else :
# 		fields_plots_color = 'ocean'
		fields_plots_color = 'seismic'
# 		fields_plots_color = 'PiYG'

	dir= os.path.dirname(folder_name+"/E/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	if spectral :
		dir= os.path.dirname(folder_name+"/E/Ex_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/E/Ex_Fourier/modulus/")
		if not os.path.exists(dir):
			os.mkdir(dir)
	
		read_and_plot.TwoD_pcolormesh_Fourier_modulus(results_folder+'/Ex_Fourier',
								N_nodes,Nx,Ny,
								'hot',
								r'$|\widehat{\mathbf{E}}_{x}|\,(m_e \Delta {\omega}^2 / e)$',
								folder_name+'/E/Ex_Fourier/modulus/Ex_F')

		dir= os.path.dirname(folder_name+"/E/Ex_Fourier/argument/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_argument(results_folder+'/Ex_Fourier_arg',
								N_nodes,Nx,Ny,
								'jet',
								r'$\mathrm{arg}\left(\widehat{\mathbf{E}}_{x}\right)\,(^o)$',
								folder_name+'/E/Ex_Fourier/argument/Ex_F_arg')

		dir= os.path.dirname(folder_name+"/E/Ey_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/E/Ey_Fourier/modulus/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_modulus(results_folder+'/Ey_Fourier',
								N_nodes,Nx,Ny,
								'hot',
								r'$|\widehat{\mathbf{E}}_{y}|\,(m_e \Delta {\omega}^2 / e)$',
								folder_name+'/E/Ey_Fourier/modulus/Ey_F')

		dir= os.path.dirname(folder_name+"/E/Ey_Fourier/argument/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_argument(results_folder+'/Ey_Fourier_arg',
								N_nodes,Nx,Ny,
								'jet',
								r'$\mathrm{arg}\left(\widehat{\mathbf{E}}_{y}\right)\,(^o)$',
								folder_name+'/E/Ey_Fourier/argument/Ey_F_arg')

		dir= os.path.dirname(folder_name+"/E/Ez_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/E/Ez_Fourier/modulus/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_modulus(results_folder+'/Ez_Fourier',
								N_nodes,Nx,Ny,
								'Spectral',#'hot',
								r'$|\widehat{\mathbf{E}}_{z}|\,(m_e \Delta {\omega}^2 / e)$',
								folder_name+'/E/Ez_Fourier/modulus/Ez_F')

		dir= os.path.dirname(folder_name+"/E/Ez_Fourier/argument/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_argument(results_folder+'/Ez_Fourier_arg',
								N_nodes,Nx,Ny,
								'Spectral',#'jet',
								r'$\mathrm{arg}\left(\widehat{\mathbf{E}}_{z}\right)\,(^o)$',
								folder_name+'/E/Ez_Fourier/argument/Ez_F_arg')

	dir= os.path.dirname(folder_name+"/E/Ez/")
	if not os.path.exists(dir):
		os.mkdir(dir)
						
	read_and_plot.TwoD_pcolormesh(results_folder+'/Ez',
							N_nodes,Nx,Ny,fields_log_scale,
							fields_plots_color,
							r'$\mathbf{E_z}\,(m_e \Delta {\omega}^2 / e)$',
							folder_name+'/E/Ez/Ez')

	dir= os.path.dirname(folder_name+"/E/Ex/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.TwoD_pcolormesh(results_folder+'/Ex',
							N_nodes,Nx,Ny,fields_log_scale,
							fields_plots_color,
							r'$\mathbf{E_x}\,(m_e \Delta {\omega}^2 / e)$',
							folder_name+'/E/Ex/Ex')

	dir= os.path.dirname(folder_name+"/E/Ey/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.TwoD_pcolormesh(results_folder+'/Ey',
							N_nodes,Nx,Ny,fields_log_scale,
							fields_plots_color,
							r'$\mathbf{E_y}\,(m_e \Delta {\omega}^2 / e)$',
							folder_name+'/E/Ey/Ey')

	dir= os.path.dirname(folder_name+"/B/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	if spectral :
		dir= os.path.dirname(folder_name+"/B/Bx_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/B/Bx_Fourier/modulus/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.TwoD_pcolormesh_Fourier_modulus(results_folder+'/Bx_Fourier',
								N_nodes,Nx,Ny,
								'Spectral',#'hot',
								r'$|\widehat{\mathbf{B}}_{x}|\,(m_e \Delta {\omega}^2 / e)$',
								folder_name+'/B/Bx_Fourier/modulus/Bx_F')

		dir= os.path.dirname(folder_name+"/B/Bx_Fourier/argument/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_argument(results_folder+'/Bx_Fourier_arg',
								N_nodes,Nx,Ny,
								'Spectral',#'jet',
								r'$\mathrm{arg}\left(\widehat{\mathbf{B}}_{x}\right)\,(^o)$',
								folder_name+'/B/Bx_Fourier/argument/Bx_F_arg')

		dir= os.path.dirname(folder_name+"/B/By_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/B/By_Fourier/modulus/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_modulus(results_folder+'/By_Fourier',
								N_nodes,Nx,Ny,
								'Spectral',#'hot',
								r'$|\widehat{\mathbf{B}}_{y}|\,(m_e \Delta {\omega}^2 / e)$',
								folder_name+'/B/By_Fourier/modulus/By_F')

		dir= os.path.dirname(folder_name+"/B/By_Fourier/argument/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_argument(results_folder+'/By_Fourier_arg',
								N_nodes,Nx,Ny,
								'Spectral',#'jet',
								r'$\mathrm{arg}\left(\widehat{\mathbf{B}}_{y}\right)\,(^o)$',
								folder_name+'/B/By_Fourier/argument/By_F_arg')

		dir= os.path.dirname(folder_name+"/B/Bz_Fourier/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+"/B/Bz_Fourier/modulus/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_modulus(results_folder+'/Bz_Fourier',
								N_nodes,Nx,Ny,
								'hot',
								r'$|\widehat{\mathbf{B}}_{z}|\,(m_e \Delta {\omega}^2 / e)$',
								folder_name+'/B/Bz_Fourier/modulus/Bz_F')
						
		dir= os.path.dirname(folder_name+"/B/Bz_Fourier/argument/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		read_and_plot.TwoD_pcolormesh_Fourier_argument(results_folder+'/Bz_Fourier_arg',
								N_nodes,Nx,Ny,
								'jet',
								r'$\mathrm{arg}\left(\widehat{\mathbf{B}}_{z}\right)\,(^o)$',
								folder_name+'/B/Bz_Fourier/argument/Bz_F_arg')

	dir= os.path.dirname(folder_name+"/B/Bx/")
	if not os.path.exists(dir):
		os.mkdir(dir)
						
	read_and_plot.TwoD_pcolormesh(results_folder+'/Bx',
							N_nodes,Nx,Ny,fields_log_scale,
							fields_plots_color,
							r'$\mathbf{B_x}\,(m_e \Delta {\omega}^2 / e)$',
							folder_name+'/B/Bx/Bx')

	dir= os.path.dirname(folder_name+"/B/By/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.TwoD_pcolormesh(results_folder+'/By',
							N_nodes,Nx,Ny,fields_log_scale,
							fields_plots_color,
							r'$\mathbf{B_y}\,(m_e \Delta {\omega}^2 / e)$',
							folder_name+'/B/By/By')

	dir= os.path.dirname(folder_name+"/B/Bz/")
	if not os.path.exists(dir):
		os.mkdir(dir)
						
	read_and_plot.TwoD_pcolormesh(results_folder+'/Bz',
							N_nodes,Nx,Ny,fields_log_scale,
							fields_plots_color,
							r'$\mathbf{B_z}\,(m_e \Delta {\omega}^2 / e)$',
							folder_name+'/B/Bz/Bz')

##########################
# Poynting #
##########################

if plot_Poynting :

	print('Poynting vector plots')
	print('----------------------------')

	if Poynting_log_scale :
		Poynting_plots_color = 'nipy_spectral'
	else :
# 		Poynting_plots_color = 'ocean'
		Poynting_plots_color = 'seismic'
# 		Poynting_plots_color = 'PiYG'

	dir= os.path.dirname(folder_name+"/Pi/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	dir= os.path.dirname(folder_name+"/Pi/Pix/")
	if not os.path.exists(dir):
		os.mkdir(dir)
						
	read_and_plot.TwoD_pcolormesh(results_folder+'/Pix',
							N_nodes,Nx,Ny,Poynting_log_scale,
							Poynting_plots_color,
							r'$\mathbf{\Pi}_x\,(m_e {\omega}^3 )$',
							folder_name+'/Pi/Pix/Pix')

	dir= os.path.dirname(folder_name+"/Pi/Piy/")
	if not os.path.exists(dir):
		os.mkdir(dir)
						
	read_and_plot.TwoD_pcolormesh(results_folder+'/Piy',
							N_nodes,Nx,Ny,Poynting_log_scale,
							Poynting_plots_color,
							r'$\mathbf{\Pi}_y\,(m_e {\omega}^3 )$',
							folder_name+'/Pi/Piy/Piy')

# 	dir= os.path.dirname(folder_name+"/Pi/Piz/")
# 	if not os.path.exists(dir):
# 		os.mkdir(dir)
# 						
# 	read_and_plot.TwoD_pcolormesh(results_folder+'/Piz',
# 							N_nodes,Nx,Ny,Poynting_log_scale,
# 							Poynting_plots_color,
# 							r'$\mathbf{\Pi}_z\,(m_e {\omega}^3 )$',
# 							folder_name+'/Pi/Piz/Piz')

####################
#  PML test  #
####################

if pml_tests_plot :

	print('PML test :')
	print('-------------')

	dir= os.path.dirname(folder_name+'/curves/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	if spectral :

		dir= os.path.dirname(folder_name+'/curves/Ez_Fourier/')
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+'/curves/Ez_Fourier/Ez_vs_x/')
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.read_file_and_plot_curve_Fourier(results_folder+'/Ez_Fourier',
								N_nodes,Nx,Ny,
								1,1,
								r'$\mathbf{E_z}$',
								folder_name+'/curves/Ez_Fourier/Ez_vs_x/Ez_vs_x')

		dir= os.path.dirname(folder_name+'/curves/Ez_Fourier/Ez_vs_y/')
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.read_file_and_plot_curve_Fourier(results_folder+'/Ez_Fourier',
								N_nodes,Nx,Ny,
								1,2,
								r'$\mathbf{E_z}$',
								folder_name+'/curves/Ez_Fourier/Ez_vs_y/Ez_vs_y')
						
		dir= os.path.dirname(folder_name+'/curves/Bx_Fourier/')
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+'/curves/Bx_Fourier/Bx_vs_x/')
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.read_file_and_plot_curve_Fourier(results_folder+'/Bx_Fourier',
								N_nodes,Nx,Ny,
								1,1,
								r'$\mathbf{B_x}$',
								folder_name+'/curves/Bx_Fourier/Bx_vs_x/Bx_vs_x')

		dir= os.path.dirname(folder_name+'/curves/Bx_Fourier/Bx_vs_y/')
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.read_file_and_plot_curve_Fourier(results_folder+'/Bx_Fourier',
								N_nodes,Nx,Ny,
								1,2,
								r'$\mathbf{B_x}$',
								folder_name+'/curves/Bx_Fourier/Bx_vs_y/Bx_vs_y')
						
		dir= os.path.dirname(folder_name+'/curves/By_Fourier/')
		if not os.path.exists(dir):
			os.mkdir(dir)

		dir= os.path.dirname(folder_name+'/curves/By_Fourier/By_vs_x/')
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.read_file_and_plot_curve_Fourier(results_folder+'/By_Fourier',
								N_nodes,Nx,Ny,
								1,1,
								r'$\mathbf{B_y}$',
								folder_name+'/curves/By_Fourier/By_vs_x/By_vs_x')

		dir= os.path.dirname(folder_name+'/curves/By_Fourier/By_vs_y/')
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.read_file_and_plot_curve_Fourier(results_folder+'/By_Fourier',
								N_nodes,Nx,Ny,
								1,2,
								r'$\mathbf{B_y}$',
								folder_name+'/curves/By_Fourier/By_vs_y/By_vs_y')

	dir= os.path.dirname(folder_name+'/curves/Ez/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	dir= os.path.dirname(folder_name+'/curves/Ez/Ez_vs_x/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.read_file_and_plot_curve(results_folder+'/Ez',
							N_nodes,Nx,Ny,
							0,1,
							r'$\mathbf{E_z}$',
							folder_name+'/curves/Ez/Ez_vs_x/Ez_vs_x')

	dir= os.path.dirname(folder_name+'/curves/Ez/Ez_vs_r/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.read_file_and_plot_curve(results_folder+'/Ez',
							N_nodes,Nx,Ny,
							1,3,
							r'$\mathbf{E_z}$',
							folder_name+'/curves/Ez/Ez_vs_r/Ez_vs_r')


	dir= os.path.dirname(folder_name+'/curves/Ez/Ez_vs_y/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.read_file_and_plot_curve(results_folder+'/Ez',
							N_nodes,Nx,Ny,
							1,2,
							r'$\mathbf{E_z}$',
							folder_name+'/curves/Ez/Ez_vs_y/Ez_vs_y')
						
	dir= os.path.dirname(folder_name+'/curves/Bx/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	dir= os.path.dirname(folder_name+'/curves/Bx/Bx_vs_r/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.read_file_and_plot_curve(results_folder+'/Bx',
							N_nodes,Nx,Ny,
							1,3,
							r'$\mathbf{B_x}$',
							folder_name+'/curves/Bx/Bx_vs_r/Bx_vs_r')

	dir= os.path.dirname(folder_name+'/curves/Bx/Bx_vs_x/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.read_file_and_plot_curve(results_folder+'/Bx',
							N_nodes,Nx,Ny,
							1,1,
							r'$\mathbf{B_x}$',
							folder_name+'/curves/Bx/Bx_vs_x/Bx_vs_x')
						
	dir= os.path.dirname(folder_name+'/curves/Bx/Bx_vs_y/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.read_file_and_plot_curve(results_folder+'/Bx',
							N_nodes,Nx,Ny,
							1,2,
							r'$\mathbf{B_x}$',
							folder_name+'/curves/Bx/Bx_vs_y/Bx_vs_y')
						
	dir= os.path.dirname(folder_name+'/curves/By/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	dir= os.path.dirname(folder_name+'/curves/By/By_vs_r/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.read_file_and_plot_curve(results_folder+'/By',
							N_nodes,Nx,Ny,
							1,3,
							r'$\mathbf{B_y}$',
							folder_name+'/curves/By/By_vs_r/By_vs_r')

	dir= os.path.dirname(folder_name+'/curves/By/By_vs_x/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.read_file_and_plot_curve(results_folder+'/By',
							N_nodes,Nx,Ny,
							1,1,
							r'$\mathbf{B_y}$',
							folder_name+'/curves/By/By_vs_x/By_vs_x')

	dir= os.path.dirname(folder_name+'/curves/By/By_vs_y/')
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.read_file_and_plot_curve(results_folder+'/By',
							N_nodes,Nx,Ny,
							1,2,
							r'$\mathbf{B_y}$',
							folder_name+'/curves/By/By_vs_y/By_vs_y')

###############
# Phase-Space #
###############

if phase_spaces_plot :

	print('Phase-space plots :')
	print('-------------------')

	dir= os.path.dirname(folder_name+"/phase_space/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	dir= os.path.dirname(folder_name+"/phase_space/ele/")
	if not os.path.exists(dir):
		os.mkdir(dir)
						
	dir= os.path.dirname(folder_name+"/phase_space/ele/px_vs_xy/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.part_plot(1,results_folder+'/ele_px_vs_xy',N_nodes,folder_name+'/phase_space/ele/px_vs_xy/ele_px_vs_xy')

	dir= os.path.dirname(folder_name+"/phase_space/ele/py_vs_xy/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.part_plot(2,results_folder+'/ele_py_vs_xy',N_nodes,folder_name+'/phase_space/ele/py_vs_xy/ele_py_vs_xy')

	dir= os.path.dirname(folder_name+"/phase_space/ele/pz_vs_xy/")
	if not os.path.exists(dir):
		os.mkdir(dir)

	read_and_plot.part_plot(3,results_folder+'/ele_pz_vs_xy',N_nodes,folder_name+'/phase_space/ele/pz_vs_xy/ele_pz_vs_xy')

	if movion :
		dir= os.path.dirname(folder_name+"/phase_space/ions/")
		if not os.path.exists(dir):
			os.mkdir(dir)
						
		dir= os.path.dirname(folder_name+"/phase_space/ions/px_vs_xy/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.part_plot(1,results_folder+'/ion_px_vs_xy',N_nodes,folder_name+'/phase_space/ions/px_vs_xy/ion_px_vs_xy')

		dir= os.path.dirname(folder_name+"/phase_space/ions/py_vs_xy/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.part_plot(2,results_folder+'/ion_py_vs_xy',N_nodes,folder_name+'/phase_space/ions/py_vs_xy/ion_py_vs_xy')

		dir= os.path.dirname(folder_name+"/phase_space/ions/pz_vs_xy/")
		if not os.path.exists(dir):
			os.mkdir(dir)

		read_and_plot.part_plot(3,results_folder+'/ion_pz_vs_xy',N_nodes,folder_name+'/phase_space/ions/pz_vs_xy/ion_pz_vs_xy')
