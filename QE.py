#!/usr/bin/env python
# -*- coding:utf-8 -*-
### options for GE job maneger
#$ -cwd
#$ -V -S /usr/local/Python-2.7.12/bin/python
#$ -N Name
#$ -e out.e
#$ -o out.o
#$ -pe smp 40
#$ -q queue_name
### options for Torque (PBS) job maneger
#$PBS -l nodes=1:ppn=12
### options for PJM job manager
#PJM -L "rscunit=ito-a"
#PJM -L "rscgrp=ito-ss"
#PJM -L "vnode=1"
#PJM -L "vnode-core=36"
#PJM -L "elapse=96:00:00"
### options for LSF job manager
#BSUB -q "queue_name"
#BSUB -n 16
#BSUB -J "Job name"
#BSUB -o std_output_file
#BSUB -m "node_name"
(T,F)=(True,False)                      #alias for bool numbers
qev=7.0                                 #qe version
#======================parallel settings===========================
(mpi_num, kthreads_num) = (12, 4)       #number of threads for MPI/openMP
ithreads_num = 0                        #number of image for ph.x
(start_q, last_q)=(0, 0)                #start and end of q to calculate in ph.x
sw_run = F                              #switch of execute DFT calculation or not
sw_mpi = T                              #switch of MPI calculation
from_poscar=False
#======================lattice parameters==========================
aa=3.218                                #lattice parameter a
ab=aa                                   #lattice parameter b
ac=5.244                                #lattice parameter c

alpha=90.                               #lattice angle alpha
beta=90.                                #lattice angle beta
gamma=120.                              #lattice angle gamma

#optional parameters
sw_celldm = F                            #use celldm(1) if ibrav==0
#cry_ax[]                                #crystal axes if you use ibrav==0
#klist=[]                                #if you choose your own klist, write it here.
#kbrav=4                                 #select build-in klist when from_poscar=True
#for supercell calculation
sw_sc=False
sc_size=[3,3,3]
imp_atom=['Eu']
imp_position=[[['Ga',0]]]
imp_pot_type=['spdn']
vac_position=[['N',0]]
#=================Crystal structure & conditions===================
prefix='GaN'                            #material name (or job name)
space=186                               #space group
atom=['Ga','N']                         #Elements Name
atomic_position=[[[1./3., 2./3., 0.],[2./3., 1./3., .5]],
                 [[1./3., 2./3., 3./8.],[2./3., 1./3., 7./8.]]]
#------------------------------------------------------------------
ibrav=4                                 #brave lattice type  
type_xc='pbe'                        #type of exchange correlation functional
pot_kind='kjpaw_psl.1.0.0'
pot_type=['dn','n']                     #psede potential name

sw_so = F                               #activate soc and noncliner calc.
sw_vdW = F                              #activate van der Waals interaction
sw_spn_pol = F
vdW_corr = 'rVV10' #'vdW-DF'            #set van der Waals type
#=================== directorys settings ==========================
sw_apw = T                              #switch of pp dir for paw( and soc) or not
outdir = './'                           #path of output directory
#mpiopt='$LSF_BINDIR/openmpi-'           #direct link of mpirun or mpiexec if you need
#================== switch of calculation =========================
sw_scf = F                              #generate input file for scf calculation
sw_dos = F                              #generate input file for dos calc.
sw_prj = F                              #generate input file for proj calc.
sw_bands = F                            #generate input file for band calc.
sw_ph = F                               #generate input file for phonon calc.
sw_dyn = F                              #generate input file f0r matdyn.x etc.
sw_epw = F                              #generate input file for epw.x 
sw_save_dir = F                         #generate save directory for epw
sw_wan = F                              #switch wannierization
sw_wan_init = F                         #generate input file for wannier90
sw_wan_init_nscf = F                    #calc. nscf cycle for wannier90
sw_post_wan =T                          #calc postwan process
sw_restart = T                          #switch restart tag or not
sw_opt = F                              #switch of optimaization
sw_md = F
#=====================pw_parameters================================
k_mesh_scf = [8,8,8]                    #k mesh for DFT calc
k_mesh_bands = 20                       #k mesh for bands calc
k_mesh_wannier = [8,8,8]                #k mesh for wannierize
ecut = 40.0                             #cut off energy of pw basis
ec_rho = 600                            #cut off energy of density
(e_conv, f_conv) = (1.0e-5, 1.0e-4)     #threshold of total energy's convergence and force's one 
(scf_conv,nscf_conv) = (1.0e-8, 1.0e-8) #threshold of convergence on scf,nscf cycles
elec_step = 200                         #threthold of scf cycle
mixing_beta = 0.5                       #mixing for scf cycle
nband = 50                              #number of bands
dgs = 0.025                             #dispersion of k-mesh
de = 0.01                               #delta E for dos
occupations='tetrahedra_opt'            #occupation setting
#occupations='smearing'
#occupations='fixed'
eband_win = [-10., 15.]                 #energy range of .ps file
#edos_win = [-30., 15.]                  #energy range of .dos file
wf_collect = T                          #collect paralleled wavefunctions or not usually T
sw_nosym = F                            #no symmetry and no inversion
opt_vol = F                             #optimize only lattice parameters
scf_mustnot_conv =F                     #noneed convergence in optimaization cycle
nspin=1                                 #spin polarized setting nonpol=1,z-axis=2,general=4
#mmom=[1,0]                              #initial spin polization
#theta_m=[90,0]                          #initial spin angle theta (tilt for z axis) use only nspin=4
#phi_m=[180,0]                           #initial spin angle phi (xy plane) use only nspin=4
#--------------------LDA+U parameters------------------------------
sw_ldaU_old = F                         #switch LDA+U (old type)
lda_U = [0., 0.]                        #Hubbard U list, length < atom
lda_J = [0., 0.]                        #Hubbard J list
#Hubb=[]                                 #Hubbard parameter setting qe 7.0<=
#------------------------HSE config--------------------------------
sw_hse=T
hse_q=[2,2,2]
#-----------------optimization cell settings-----------------------
opt_step = 100                          #number of optimization step
press=0.0                               #pressure (Kbar)
p_conv=.5e0                             #conv threshold
#----------------------MD settings---------------------------------
md_step=500                             #number of MD step
dt = 20                                 #time step 20a.u.~1fs
sw_vc = F                               #variable lattice param. or not
dynamics = 'verlet'
md_temp = 300                           #temperature for MD
ion_temp = 'andersen'
#======================ph_parameters===============================
q_mesh_dyn = [4, 4, 4]                  #q mesh for phonon DFPT calc
q_mesh_bands = 20                       #q mesh for phonon band calc
q_mesh_dos = 8                          #q mesh for phonon dos calc
ph_conv = 1.0e-14                       #threshold of energy's convergence for phonon
amix = 0.7                              #mixing rate for ph scf default=0.7
maxiter_ph = 100                        #max iteration default=100
pband_win = [0, 200]                    #energy range of .ps file
sw_ep = F                               #swich to calc. e-p interaction or not we cannot obtain e-p int. in nonclin.
mpol = F                                #if material is local polarized mpol=T
qnum = 10                               #number of q smearing
sw_qshift = T                           #if calc. lambda with opt_tetra, make this T
sw_gen_a2f = F                          #switch to generage a2f.dat files
#===================Wannier_parameters=============================
# the projections name which we can use are s,p,d,f, and sp,sp2,sp3, sp3d,sp3d2
nwann = 12                              #number of wannier basis
dis_win = [-10.00, 15.00]               #max(min)_window, range of sub space energy
frz_win = [-10.0, 6.0]                  #froz_window dp22
projection = [('Ga','p'),('N','p')]     #projections, initial funcution of wannier
sw_fs_plot = F                          #plot Fermi surface
fermi_mesh = 100                        #mesh of k-points in bxsf file
unk = F                                 #Bloch(Wannier)_func
uwrite = F                              #output unitary matrix for bloch to wannier
#=====================postwan_parameters===========================
#---------------------Bolz_wann-----------------------------------
boltz_kmesh=40                          #k-mesh size of boltzwann
btau=1                                  #relaxation time of boltwann
mu_range=[0.,0.]                        #mu range of boltzwann
bmu_step=1                              #mu step of boltzwann
temp_range=[50,1000]                    #temp range of boltzwann
bt_step=50                              #temp step size of boltzwann
#=================EPW parameters===================================
epw_k_mesh=[16,16,16]                   #interpolate k-mesh size for epw
epw_q_mesh=[8,8,8]                      #interpolate q-mesh size for epw
#=============================modules==============================
import numpy as np
import os, datetime, subprocess, argparse
#=====================global_lambda_expression=====================
TorF=lambda x:'.True.' if x else '.False.'
w_conv=lambda a:'1.0E%d'%int(np.log10(a))
#====================== physical parameters =======================
sig_fig = 9                             #significant figure after decimal fraction
bohr=round(0.52917721067, sig_fig)      #Bohr Radius
ibohr=1.0/bohr                          #inverse of Bohr Radius
#========================= atomic mass ============================
#       1             2             13            14          15          16           17            18
mass={'H':1.00794,                                                                                  'He':4.002602,
      'Li':6.941,    'Be':9.012182,'B':10.811,   'C':12.0107,'N':14.0067,'O':15.9994,'F':18.9984032,'Ne':20.1797,
      'Na':22.98977, 'Mg':24.305,  'Al':26.981538,'Si':28.0855,'P':30.973761,'S':32.065,'Cl':35.453,'Ar':39.948,
      'K':39.0983,   'Ca':40.078,
                               'Sc':44.95591,'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049, #3d
                               'Fe':55.845,'Co':58.9332,'Ni':58.6934,'Cu':63.546,'Zn':65.409,     #3d
                                   'Ga':69.723,  'Ge':72.64, 'As':74.9216,'Se':78.96, 'Br':79.904,  'Kr':83.798,
      'Rb':85.4678,  'Sr':87.62,
                               'Y':88.90585,'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':97.907216,  #4d
                               'Ru':101.07,'Rh':102.9055,'Pd':106.42,'Ag':107.8682,'Cd':112.411,  #4d
                                   'In':114.818, 'Sn':118.71,'Sb':121.76,'Te':127.6,  'I':126.90447,'Xe':131.293,
      'Cs':132.90545,'Ba':137.327,
                               'La':138.9055, #lanthanide
             'Ce':140.116,'Pr':140.90765,'Nd':144.24,'Pm':144.912744,'Sm':150.36,'Eu':151.964,'Gd':157.25,
             'Tb':158.92534,'Dy':162.5,'Ho':164.93032,'Er':167.259,'Tm':168.93421,'Yb':173.04,'Lu':174.967,
                               'Ln':138.9055,'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,   #5d
                               'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,'Hg':200.59,  #5d
                                   'Tl':204.3833,'Pb':207.2,'Bi':208.98038,
             'Th':232.0381,'Pa':231.03588,'U':238.02891} #actinide
#==========================global_variables========================
axis=[aa,ab,ac]                         #lattice parameters a,b,c
deg=[alpha,beta,gamma]                  #lattice parameters alpha,beta,gamma
axis=np.array(axis)
fildvscf="'dvscf'"
fildyn='%s.dyn'%prefix
flfrc='%s.fc'%prefix
recover=TorF(sw_restart)
dxml='' #'.xml' if sw_so else '' #some old version add .xml when with soc phonon calculation
#==============serch several variables and initialize if can't find
parser=argparse.ArgumentParser(prog='QE.py',description='input file generator for Quantum Espresso')
parser.add_argument("-s","-scf",help='generate scf file',action='store_true')
parser.add_argument("-d","-dos",help='generate dos calc files',action='store_true')
parser.add_argument("-b","-band",help='generate band calc files',action='store_true')
parser.add_argument("-p","-ph","-phonon",help='generate phonon calc files',action='store_true')
parser.add_argument("-w","-wan","-wannier",help='generate wannier calc files',action='store_true')
parser.add_argument("-o","-opt","-optimize",help='generate scf file for optimization',action='store_true')
parser.add_argument("-e","-epw","-el_phonon",help='generate epw calc files',action='store_true')
parser.add_argument("-m","-md",help='generate scf file for MD',action='store_true')
args=parser.parse_args()
if args.s:
    sw_scf=True
if args.d:
    sw_dos=True
if args.b:
    sw_bands=True
if args.p:
    sw_ph=True
    sw_dyn=True
if args.w:
    sw_wan=True
    sw_wan_init = True
    sw_wan_init_nscf = True
if args.o:
    sw_opt=True
if args.e:
    sw_epw=True
if args.m:
    sw_md=True
try:
    mpiopt
except NameError:
    mpiopt=''
try:
    edos_win
except NameError:
    edos_win=eband_win

if sw_fs_plot:
    try: #detect fermi_mesh
        fermi_mesh
    except NameError:
        fermi_mesh=100

#==========================functions===============================
def get_bravs(space,ibrav):
    if ibrav==0:
        num_brav=0
    else:
        if isinstance(space,int):
            if space in {23,24,44,45,46,71,72,73,74,79,80,82,87,88,
                         97,98,107,108,109,110,119,120,121,122,139,
                         140,141,142,197,199,204,206,211,214,217,220,229,230}: #I (Body Center)
                if space >195: #BCC
                    num_brav=3
                elif space>75: #BCT
                    num_brav=7
                else: #BCO
                    num_brav=11
            elif space in {22,42,43,69,70,196,202,203,209,210,216,219,225,226,227,228}: #F (Face Center)
                if space >195: #FCC
                    num_brav=2
                else: #FCO
                    num_brav=10
            elif space in {146,148,155,160,161,166,167}: #R (trigonal(rhombohedral))
                num_brav=5
            elif space in {38,39,40,41,5,8,9,12,15,20,21,35,36,37,63,64,65,66,67,68}: #ABC (Base Center)
                if space >16: #Ortho
                    num_brav=9
                else: #monocli
                    if deg[2]!=90.:
                        num_brav=13
                    elif deg[1]!=90.:
                        num_brav=-13
            else: #P (Simple)
                if space >194: #SC
                    num_brav=1
                elif space>142: #Trigonal or Hexagonal
                    num_brav=4
                elif space>75: #ST
                    num_brav=6
                elif space>15: #SO
                    num_brav=8
                elif space>2: #monocli
                    if deg[2]!=90.:
                        num_brav=12
                    elif deg[1]!=90.:
                        num_brav=-12
                else:
                    num_brav=14
        elif isinstance(space,str):
            if 'I' in space: #Body Center
                if '3' in space: #Cube
                    num_brav=3
                elif '4' in space: #Tetra
                    num_brav=7
                else: #Ortho
                    num_brav=11
            elif 'F' in space: #Face Center
                if '3' in space: #Cube
                    num_brav=2
                else: #Ortho
                    num_brav=10
            elif ('R' in space): #Trigonal
                num_brav=5
            elif ('A' in space): #Base Center
                pass
            elif ('C' in space): #Base Center
                num_brav=9
            else: #Simple
                if '6' in space: #Hexagonal
                    num_brav=4
                elif '3' in space:
                    if deg[2]==90: #Cube
                        num_brav=1
                    else: #Trigonal
                        num_brav=4
                elif '4' in space: #Tetra
                    num_brav=6
                elif '2' in space or 'm' in space:
                    ck=(space.count('2')+space.count('m')
                        +space.count('n')+space.count('c')
                        +space.count('a')+space.count('b'))
                    if ck==3: #Ortho
                        num_brav=8
                    elif ck==2:
                        num_brav=8
                    else:
                        num_brav=13
                else:
                    num_brav=14
    ibrav=num_brav
    return ibrav,num_brav

def gen_SC_positions(sc_size,positions):
    sc_positions=[]
    for pos in positions:
        tmp=[]
        for i in range(sc_size[2]):
            z_slide=1./sc_size[2]
            for p in pos:
                for j in range(sc_size[1]):
                    y_slide=1./sc_size[1]
                    for k in range(sc_size[0]):
                        x_slide=1./sc_size[0]
                        slide=np.array([k*x_slide,j*y_slide,i*z_slide])
                        tmp.append(list(p*np.array([x_slide,y_slide,z_slide])+slide))
        sc_positions.append(tmp)
    try:
        vac_position
        del_num=[0]*len(atom)
        for pos in vac_position:
            for i,at in enumerate(atom):
                if at==pos[0]:
                    del sc_positions[i][pos[1]-del_num[i]]
                    del_num[i]+=1
    except NameError:
        print('no vacancy',flush=True)
    try:
        imp_atom
        imp_position
        if len(imp_atom)!=0 and len(imp_atom)==len(imp_position):
            new_imp_pos=[]
            for pos in imp_position:
                tmp=[]
                for ps in pos:
                    for i,at in enumerate(atom):
                        if at==ps[0]:
                            tmp.append(sc_positions[i].pop(ps[1]))
                new_imp_pos.append(tmp)
        sc_positions=new_imp_pos+sc_positions
    except NameError:
        print('no impurity',flush=True)
    return sc_positions

def read_poscar(fname='POSCAR'):
    f=open(fname,'r')
    data=f.readlines()
    a=float(data[1])
    cry_ax0=np.array([[float(d) for d in dd.strip().split()] for dd in data[2:5]])
    axis=a*np.sqrt((abs(cry_ax0)**2).sum(axis=1))
    cry_ax=a*cry_ax0/axis
    atom=data[5].strip().split()
    na=[int(a) for a in data[6].strip().split()]
    atomic_position=[]
    cons=8
    for i in na:
        tmp=[[float(d) for d in dl.strip().split()] for dl in data[cons:cons+i]]
        atomic_position.append(tmp)
        cons+=i
    return(axis,cry_ax,atom,atomic_position)

def generate_klist(num_brav):
    if num_brav==1: #simple cube
        k_list=[['R',[.5,.5,0.]],['G',[0.,0.,0.]],
                ['X',[.5,0.,0.]],['M',[.5,.5,0.]],
                ['G',[0.,0.,0.]]]
    elif num_brav==2: #fcc
        k_list=[['G',[0.,0.,0.]],['X',[0.5,0.,0.5]],
                ['G',[1.0,0.,0.]],['L',[0.5,0.5,0.5]],
                ['W',[0.5,0.25,0.75]],['G',[0.,0.,0.]]]
    elif num_brav==3: #bcc
        k_list=[['G',[0.,0.,0.]],['H',[0.5,0.5,0.5]],['N',[0.,0.5,0.]],
                ['P',[0.25,0.75,-0.25]],['G',[0.,0.,0.]],['N',[0.,0.5,0.]]]
    elif num_brav==4: #hexagonal
        k_list=[['G',[0.,0.,0.]],['K',[2./3,-1./3,0.]],
                ['M',[.5,0.,0.]],['G',[0.,0.,0.]],['Z',[0.,0.,.5]]]
    elif num_brav==5: #trigonal
        k_list=[['G',[0.,0.,0.]],['K',[2./3,-1./3,0.]],
                ['M',[.5,0.,0.]],['G',[0.,0.,0.]],['Z',[0.,0.,.5]]]
    elif num_brav==6: #simple tetra
        k_list=[['G',[0.,0.,0.]],['X',[.5,0.,0.]],
                ['M',[.5,.5,0.]],['G',[0.,0.,0.]],
                ['Z',[0.,0.,.5]]]
    elif num_brav==7: #bct
        k_list=[['G',[0.,0.,0.]],['X',[0.0,0.5,-0.5]],
                ['P',[0.25,0.75,-0.25]],['N',[0.,0.5,0.]],
                ['G',[0.,0.,0.]],['Z',[0.5,0.5,0.5]]]
    elif num_brav==8: #simple orthogonal
        k_list=[['G',[0.,0.,0.]],['X',[0.5,0.0,0.]],
                ['M',[0.5,0.5,0.]],['Y',[0.,0.5,0.]],
                ['G',[0.,0.,0.]],['M',[0.5,0.5,0.]]]
    elif num_brav==9: #base centerd ortho
        k_list=[['G',[0.,0.,0.]],['X',[.5,0.,0.]],['M',[.5,.5,0.]],
                ['Y',[0.,.5,0.]],['G',[0.,0.,0.]],['M',[.5,.5,0.]]]
    elif num_brav==10: #fco
        k_list=[['G',[0.,0.,0.]],['X',[.5,0.,0.]],['M',[.5,.5,0.]],
                ['Y',[0.,.5,0.]],['G',[0.,0.,0.]],['M',[.5,.5,0.]]]
    elif num_brav==11: #bco
        k_list=[['G',[0.,0.,0.]],['X',[.5,0.,0.]],['M',[.5,.5,0.]],
                ['Y',[0.,.5,0.]],['G',[0.,0.,0.]],['M',[.5,.5,0.]]]
    elif num_brav==12: #monocli
        k_list=[['G',[0.,0.,0.]],['X',[.5,0.,0.]],['M',[.5,.5,0.]],
                ['Y',[0.,.5,0.]],['G',[0.,0.,0.]],['M',[.5,.5,0.]]]
    elif num_brav==13: #monocli
        k_list=[['G',[0.,0.,0.]],['X',[.5,0.,0.]],['M',[.5,.5,0.]],
                ['Y',[0.,.5,0.]],['G',[0.,0.,0.]],['M',[.5,.5,0.]]]
    else: #monocli
        k_list=[['G',[0.,0.,0.]],['X',[.5,0.,0.]],['M',[.5,.5,0.]],
                ['Y',[0.,.5,0.]],['G',[0.,0.,0.]],['Z',[0.,0.,.5]]]
    return k_list

def write_file(name,stream):
    f=open(name,'w')
    f.write(stream)
    f.close()

def check_type(obj,typ):
    if not isinstance(obj,typ):
        print('type err')
        exit()
    else:
        pass

def date():
    from locale import setlocale, LC_ALL
    d=datetime.datetime.today()
    setlocale(LC_ALL,'')
    print(d.strftime('%Y年  %B %d日 %A %H:%M:%S JST' if os.environ['LANG'].find('ja')!=-1 
                     else '%a %b %d %X JST %Y'),flush=True)

def os_and_print(command):
    print(command,flush=True)
    info=subprocess.run(command,shell=True)
    if info.returncode!=0:
        print('something error',flush=True)
        exit()

def get_ef(name,ext):
    fname='%s.%s.out'%(name,ext)
    if os.path.exists(fname):
        if occupations=='fixed':
            ckstrings='lowest unoccupied level'
            spst='(ev):'
            spst2=''
        else:
            ckstrings='Fermi'
            spst='is'
            spst2='ev'
        for f in open(fname,'r'):
            if f.find(ckstrings)!=-1:
                item=f.split(spst)
                it1=item[1].split(spst2)
                return float(it1[0])
        else:
            print('input error from %s.%s.out\n'%(name,ext)
                  +'return ef = 0\n',flush=True)
            return 0.
    else:
        print('can not find %s\n return ef = 0\n'%fname,flush=True)
        return 0.

def make_fstring_obj(obj_name,var_list,val_dic,sw_form):
    if sw_form=='pw':
        form='%28s = '
    elif sw_form=='ph':
        form='%20s = '
    elif sw_form=='pw2wan':
        form='%15s = '
    else:
        form='%12s = '
    fstring='&%s\n'%obj_name
    for vc in var_list:
        fstring+=form%vc+str(val_dic[vc])+'\n'
    fstring+='/\n'
    return fstring

def atom_position(atom,atomic_position):
    mat=get_cr_mat(num_brav,F)
    mat=np.linalg.inv(mat).T
    aposition=[[list(mat.dot(np.array(ap))) for ap in app] for app in atomic_position]
    tmp='%2s  %12.9f %12.9f %12.9f'+(' 0 0 0\n' if(sw_opt and opt_vol) else '\n')
    atom_string=''
    for i,at in enumerate(atom):
        for ap in aposition[i]:
            atom_string+=tmp%tuple([at]+ap)

    return atom_string

def get_cr_mat(num_brav,sw=T):
    if num_brav in {1,6,8}: #Simple
        mat=np.identity(3)
    elif num_brav in {2,10}: #Face center
        mat=np.array([[-.5, 0., .5],
                      [ 0., .5, .5],
                      [-.5, .5, 0.]])
    elif num_brav in {3,7,11}: #Body center
        mat=np.array([[ .5, -.5, .5],
                      [ .5,  .5, .5],
                      [-.5, -.5, .5]])
    elif num_brav==4: #hexa
        if sw:
            mat=np.array([[ 1., 0.            , 0.],
                          [-.5, .5*np.sqrt(3.), 0.],
                          [ 0., 0.            , 1.]])
        else:
            mat=np.identity(3)
    elif num_brav==5: #trigonal
        cg=np.cos(np.pi*deg[2]/180)
        tx=np.sqrt((1-cg)*0.5)
        ty=np.sqrt((1-cg)/6)
        tz=np.sqrt((1+2*cg)/3)
        if sw:
            mat=np.array([[ tx,  -ty, tz],
                          [ 0., 2*ty, tz],
                          [-tx,  -ty, tz]])
        else:
            mat=np.identity(3)
    elif num_brav==9: #Orthorhombic Base center
        mat=np.array([[ .5, .5, 0.],
                      [-.5, .5, 0.],
                      [ 0., 0., 1.]])
    elif num_brav==12: #Monoclinic
        pass
    elif num_brav==13: #Monoclinic
        pass
    elif num_brav==14: #Monoclinic
        if sw:
            phase=np.pi*deg[0]/180.
            phase2=np.pi*deg[1]/180.
            phase3=np.pi*deg[2]/180.
            ax1=axis[1]/axis[0]
            ax2=axis[2]/axis[0]
            r1=ax1*np.cos(phase3)
            r2=ax1*np.sin(phase3)
            r3=ax2*np.cos(phase2)
            r4=ax2*(np.cos(phase)-np.cos(phase2)*np.cos(phase3))/np.sin(phase3)
            r5=ax2*np.sqrt(1+2*np.cos(phase)*np.cos(phase2)*np.cos(phase3)
                           -(np.cos(phase)**2+np.cos(phase2)**2+np.cos(phase3)**2))/np.sin(phase3)
            mat=np.array([[ 1., 0.,  0.],
                          [ r1, r2,  0.],
                          [ r3, r4, r5]])
        else:
            mat=np.identity(3)
    else:
        try:
            cry_ax
            if sw:
                mat=np.array(cry_ax)
            else:
                mat=np.identity(3)
        except NameError:
            print('please set crystal axes matrix cry_ax')
            exit()
    return mat

def cell_parameter_stream(axis,deg):
    mat=get_cr_mat(num_brav)
    if ibrav==0 and sw_celldm:
        a_vec=list(mat)
    else:
        mat_ax=np.identity(3)*axis*ibohr
        a_vec=list(mat.dot(mat_ax))
    cell_string=''
    for aa in a_vec:
        cell_string=cell_string+'  %12.9f  %12.9f  %12.9f\n'%tuple(aa)

    return cell_string

def k_line_stream(k_num,k_list):
    k_string='%d\n'%(k_num*(len(k_list)-1)+1)
    dk=1./k_num
    w=1.
    wt=' %f'%w if True else ''
    for kb,ka in zip(k_list,k_list[1:]):
        ks=kb[1]
        kp=[(kaa-kbb)*dk for kaa,kbb in zip(ka[1],kb[1])]
        for i in range(k_num):
            k_string+='%11.8f %11.8f %11.8f'%(ks[0]+kp[0]*i,ks[1]+kp[1]*i,ks[2]+kp[2]*i)+wt+'\n'
    k_string+='%11.8f %11.8f %11.8f'%tuple(k_list[-1][1])+wt+'\n'
    return k_string

def k_cube_stream(k_num,w_sw,sw_wan):
    wfunc=lambda x,y,z:'\n' if x else '  %10.8f\n'%(1./z if y else 1.)
    if not isinstance(w_sw,bool):
        w_sw=False
    if isinstance(k_num,int):
        dk=1./k_num
        wst=wfunc(sw_wan,w_sw,(k_num**3))
        k_string='' if sw_wan else '%d\n'%k_num**3
        for i in range(k_num):
            for j in range(k_num):
                for k in range(k_num):
                    k_string+='  %10.8f  %10.8f  %10.8f'%(dk*i,dk*j,dk*k)+wst
    elif isinstance(k_num,list):
        if len(k_num)==1:
            dk=[1./knum[0]]*3
            kn=k_num*3
        if len(k_num)==2:
            dk=[1./kp for kp in k_num]
            dk=[dk[0]]+dk
            kn=[k_num[0]]+k_num
        elif len(k_num)==3:
            dk=[1./kp for kp in k_num]
            kn=k_num
        else:
            print('k dimension <= 3')
            exit()
        w0=1.
        if w_sw:
            for wi in kn:
                w0=w0*wi
        wst=wfunc(sw_wan,w_sw,w0)
        k_string='' if sw_wan else '%d\n'%w0
        for i in range(kn[0]):
            for j in range(kn[1]):
                for k in range(kn[2]):
                    k_string+='  %10.8f  %10.8f  %10.8f'%(dk[0]*i,dk[1]*j,dk[2]*k)+wst
    else:
        print('Please input list or int into k_point ')
        exit()

    return k_string

#---------------------input file generators-------------------------------
def make_pw_cp_in(calc,kconfig,restart="'from_scratch'"):
    def atomic_parameters_stream(atom,atomic_position,UPF):
        atom_string='\nATOMIC_SPECIES\n'
        for at,up in zip(atom,UPF):
            atom_string+=' %-2s %11.7f  %s\n'%(at,mass[at[:2]],up)
        atom_string+='\nATOMIC_POSITIONS crystal\n'
        atom_string+=atom_position(atom,atomic_position)
        atom_string+='\n'
        return atom_string

    (fext,convthr) =(('nscf' if calc=='nscf' else 'bands',nscf_conv) 
                     if calc in {'nscf','bands'} else ('md',scf_conv) if calc in {'md','vc-md'}
                     else('scf',scf_conv))
    occup=("'smearing'" if kconfig else "'%s'"%occupations)

    fname='%s.%s'%(prefix,fext)
    fstream=''
    var_control=['title','calculation','restart_mode','outdir','pseudo_dir',
                 'prefix','etot_conv_thr','forc_conv_thr','wf_collect']
    if not (sw_ldaU_old and (sum(lda_J)!=0 or sw_so)):
        var_control=var_control+['tstress','tprnfor']
    val_control={'title':"'%s'"%prefix,'calculation':"'%s'"%calc,'restart_mode':restart,'outdir':"'%s'"%outdir,
                 'pseudo_dir':"'%s'"%pseude_dir,'prefix':"'%s'"%prefix,'etot_conv_thr':w_conv(e_conv),
                 'forc_conv_thr':w_conv(f_conv),'tstress':'.True.','tprnfor':'.True.',
                 'wf_collect':TorF(wf_collect),'verbosity':'high'}
    if calc in {'relax','md','vc-relax','vc-md'}:
        nstep=opt_step if calc in {'relax','vc-relax'} else md_step
        var_control+=['nstep']
        val_control.update({'nstep':nstep})
        if calc in {'md','vc-md'}:
            var_control+=['dt']
            val_control.update({'dt':dt})
    fs_control=make_fstring_obj('control',var_control,val_control,'pw')
    fstream+=fs_control

    var_system=['ibrav','nat','ntyp','occupations','smearing','degauss','nbnd','ecutwfc']
    val_system={'ibrav':ibrav,'nat':sum(len(a) for a in atomic_position),
                'ntyp':len(atom),'occupations':occup,'smearing':"'marzari-vanderbilt'",
                'degauss':0.025,'la2f':'.True.','nbnd':nband,'ecutwfc':ecut,'nosym':'.True.',
                'noinv':'.True.','noncolin':'.True.','lspinorb':'.True.','nspin':nspin}
    if sw_ep:
        var_system+=['la2f']
    if sw_nosym:
        var_system+=['nosym','noinv']
    try:
        ec_rho
        var_system+=['ecutrho']
        val_system.update({'ecutrho':ec_rho})
    except NameError:
        pass
    if sw_so:
        var_system+=['noncolin','lspinorb']
    if sw_vdW:
        if vdW_corr=='vdW-DF':
            if sw_so:
                print('We cannot use vdW-DF with noncllinear spin',flush=True)
            else:
                var_system+=['input_dft']
                val_system.update({'input_dft':"'vdW-DF'"})
        elif vdW_corr=='rVV10':
            if sw_so:
                print('We cannot use rVV10 with noncllinear spin',flush=True)
            else:
                var_system+=['input_dft']
                val_system.update({'input_dft':"'rVV10'"})
        else:
            var_system+=['vdw_corr']
            val_system.update({'vdw_corr':"'%s'"%vdW_corr})
    else:
        if sw_hse and calc=='scf':
            var_system+=['input_dft','nqx1','nqx2','nqx3','x_gamma_extrapolation','exxdiv_treatment']
            val_system.update({'input_dft':"'hse'",'nqx1':hse_q[0],'nqx2':hse_q[1],'nqx3':hse_q[2],
                               'x_gamma_extrapolation':'.True.','exxdiv_treatment':'gygi-baldereschi'})
    if sw_spn_pol:
        var_system+=['nspin']
        if len(mmom)!=0:
            for i,m in enumerate(mmom):
                var_system+=['starting_magnetization(%d)'%(i+1)]
                val_system.update({'starting_magnetization(%d)'%(i+1):m})
            if nspin==4:
                if len(theta_m)!=0:
                    for i,m in enumerate(theta_m):
                        var_system+=['angle1(%d)'%(i+1)]
                        val_system.update({'angle1(%d)'%(i+1):m})
                if len(phi_m)!=0:
                    for i,m in enumerate(phi_m):
                        var_system+=['angle2(%d)'%(i+1)]
                        val_system.update({'angle2(%d)'%(i+1):m})
    if sw_ldaU_old:
        var_system+=['lda_plus_u']
        val_system.update({'lda_plus_u':'.True.'})
        if qev<7.0:
            if sw_so or sum(lda_J)!=0:
                var_system+=['lda_plus_u_kind']
                val_system.update({'lda_plus_u_kind':1})
            for i,U in enumerate(lda_U[:len(atom)]):
                tmp='Hubbard_U(%d)'%(i+1)
                var_system+=[tmp]
                val_system.update({tmp:U})
            if sum(lda_J)!=0:
                for i,J in enumerate(lda_J[:len(atom)]):
                    tmp='Hubbard_J(1,%d)'%(i+1)
                    var_system+=[tmp]
                    val_system.update({tmp:J})
        else:
            var_system+=['lda_plus_u_kind']
            val_system.update({'lda_plus_u_kind':0})
            for i,U in enumerate(lda_U[:len(atom)]):
                tmp='Hubbard_alpha(%d)'%(i+1)
                var_system+=[tmp]
                val_system.update({tmp:U})
            if sum(lda_J)!=0:
                for i,J in enumerate(lda_J[:len(atom)]):
                    tmp='Hubbard_beta(%d)'%(i+1)
                    var_system+=[tmp]
                    val_system.update({tmp:J})
    if ibrav!=0:
        """
        setting cell parameters
        ibrav = 1~3: cube (1: P, 2: F, 3: I)
              = 4: Hexagonal and Trigonal P
              = 5,-5: Trigonal R 3fold axis 5: c or -5: <1,1,1>
              = 6,7: Tetragnal (6: P, 7: I)
              = 8~11: Orthorhombic (8: P, 9,-9: B, 10: F, 11: I)
              = 12~13: Monoclinic (12,-12: P, 13: B)
              = 14: Triclinic
        """
        var_system+=['celldm(1)'] 
        val_system.update({'celldm(1)':round(axis[0]*ibohr,sig_fig)})
        if ibrav in {4,6,7,8,9,10,11,12,13,14,-12,-13}:
            if ibrav in {8,9,10,11,12,13,14,-12}:
                var_system+=['celldm(2)']
                val_system.update({'celldm(2)':round(axis[1]/axis[0],sig_fig)})
            var_system+=['celldm(3)']
            val_system.update({'celldm(3)':round(axis[2]/axis[0],sig_fig)})
            if ibrav in {12,13}:
                var_system+=['celldm(4)']
                val_system.update({'celldm(4)':round(np.cos(np.pi*deg[2]/180.),sig_fig)})
            elif ibrav in {-12,-13}:
                var_system+=['celldm(5)']
                val_system.update({'celldm(5)':round(np.cos(np.pi*deg[1]/180.),sig_fig)})
            elif ibrav==14:
                var_system+=['celldm(4)','celldm(5)','celldm(6)']
                val_system.update({'celldm(4)':round(np.cos(np.pi*deg[0]/180.),sig_fig),
                                   'celldm(5)':round(np.cos(np.pi*deg[1]/180.),sig_fig),
                                   'celldm(6)':round(np.cos(np.pi*deg[2]/180.),sig_fig)})
        elif ibrav in (5,-5):
            var_system+=['celldm(4)']
            val_system.update({'celldm(4)':round(np.cos(np.pi*deg[2]/180.),sig_fig)})
    else:
        if sw_celldm:
            var_system+=['celldm(1)']
            val_system.update({'celldm(1)':round(axis[0]*ibohr,sig_fig)})
    fs_system=make_fstring_obj('system',var_system,val_system,'pw')
    fstream+=fs_system

    var_electrons=['diagonalization','conv_thr']
    val_electrons={'diagonalization':"'david'",'conv_thr':w_conv(convthr),'scf_must_converge':'.False.'}
    try:
        elec_step
        var_electrons+=['electron_maxstep']
        val_electrons.update({'electron_maxstep':elec_step})
    except NameError:
        pass
    try:
        mixing_beta
        var_electrons+=['mixing_beta']
        val_electrons.update({'mixing_beta':mixing_beta})
    except NameError:
        pass
    if sw_opt:
        if scf_mustnot_conv:
            var_electrons+=['scf_must_converge']
    fs_electrons=make_fstring_obj('electrons',var_electrons,val_electrons,'pw')
    fstream+=fs_electrons

    if calc in {'relax','md','vc-relax','vc-md'}:
        var_ions=[]
        val_ions={}
        if calc in {'md','vc-md'}:
            var_ions+=['ion_dynamics','ion_temperature','tempw']
            val_ions.update({'ion_dynamics':dynamics,'ion_temperature':ion_temp,'tempw':md_temp})
        fs_ions=make_fstring_obj('ions',var_ions,val_ions,'pw')
        fstream+=fs_ions

    if calc in {'vc-relax','vc-md'}:
        var_cell=[]
        val_cell={}
        if press!=0:
            var_cell+=['press','press_conv_thr']
            val_cell.update({'press':press,'press_conv_thr':p_conv})
        fs_cell=make_fstring_obj('cell',var_cell,val_cell,'pw')
        fstream+=fs_cell

    fs_atomparam=atomic_parameters_stream(atom,atomic_position,UPF)
    fstream+=fs_atomparam
    if ibrav==0:
        fs_cellparam='CELL_PARAMETERS\n'+cell_parameter_stream(axis,deg)+'\n'
        fstream+=fs_cellparam
    fstream+='K_POINTS %s\n'%('CRYSTAL' if kconfig else 'AUTOMATIC')
    if kconfig:
        if calc=='nscf':
            fs_kl_param=k_cube_stream(k_mesh_wannier,T,F)
        else:
            fs_kl_param=k_line_stream(k_mesh_bands,k_list)
        fstream+=fs_kl_param
    else:
        if isinstance(k_mesh_scf,list):
            if len(k_mesh_scf)==3:
                k_mesh=k_mesh_scf
            elif len(k_mesh_scf)==2:
                k_mesh=[k_mesh_scf[0]]*2+[k_mesh_scf[1]]
            else:
                k_mesh=[k_mesh_scf[0]]*3
        elif isinstance(k_mesh_scf,int):
            k_mesh=[k_mesh_scf]*3
        else:
            print('k_mesh_scf is list(datatype=int) or int only')
            exit()
        fstream+='%d %d %d %d %d %d\n'%tuple(k_mesh+[0]*3)
    try:
        Hubb
        if len(Hubb)!=0:
            fshubb='\nHUBBARD atomic\n'
            for hub in Hubb:
                fshubb+='  %s %s %5.2f\n'%tuple(hub)
            fstream+=fshubb
    except NameError:
        pass
    write_file(fname,fstream)

def make_dos_in():
    fname='%s.dos_in'%prefix
    fildos='%s.dos'%prefix
    ef=get_ef(prefix,'scf')
    var_dos=['prefix','outdir','Emin','Emax','DeltaE','ngauss','degauss','fildos']
    val_dos={'prefix':"'%s'"%prefix,'outdir':"'%s'"%outdir,'fildos':"'%s'"%fildos,
               'Emin':edos_win[0]+ef,'Emax':edos_win[1]+ef,'DeltaE':de,'ngauss':1,'degauss':2*dgs}
    fstream=''
    fs_dos=make_fstring_obj('dos',var_dos,val_dos,'dos')
    fstream+=fs_dos
    write_file(fname,fstream)

def make_prjwfc_in():
    fname='%s.prjwfc'%prefix
    filpdos='%s.pdos'%prefix
    var_prj=['prefix','outdir','Emin','Emax','DeltaE','ngauss','degauss','filpdos','pawproj']
    val_prj={'prefix':prefix,'outdir':"'%s'"%outdir,'filpdos':"'%s'"%filpdos,
             'Emin':edos_win[0],'Emax':edos_win[1],'DeltaE':de,'ngauss':1,'degauss':dgs,
             'pawproj':'.Ture.'}
    fstream=''
    fs_prj=make_fstring_obj('projwfc',var_prj,val_prj,'prjwfc')
    fstream+=fs_prj
    write_file(fname,fstream)

def make_bands_in():
    fname='%s.bands_in'%prefix
    fstream=''
    filband='%s_bands.dat'%prefix
    var_bands=['prefix','outdir','filband']
    val_bands={'prefix':prefix,'outdir':"'%s'"%outdir,'filband':filband}
    fs_bands=make_fstring_obj('bands',var_bands,val_bands,'bands')
    fstream=fstream+fs_bands
    write_file(fname,fstream)

def make_pw2wan_in():
    fname='%s.pw2wan'%prefix
    fstream='\n'
    var_bands=['outdir','prefix','seedname','spin_component','write_unk']
    val_bands={'outdir':"'./'",'prefix':"'%s'"%prefix,'seedname':"'%s'"%prefix,
               'spin_component':"'none'",'write_unk': TorF(unk)}
    fs_bands=make_fstring_obj('inputpp',var_bands,val_bands,'pw2wan')
    fstream+=fs_bands
    write_file(fname,fstream)

def make_win():
    def win_strings(values,sname):
        def format_val(a):
            if isinstance(a,float):
                return '%7.4f'%a
            else:
                return str(a)
        if(sname=='plot'):
            form='%-25s = '
        else:
            form='%-12s = '
        strings=''
        for vs in values:
            strings+=form%vs[0]+format_val(vs[1])+'\n'
        strings+='\n'
        return strings
    def get_k_point_path():
        p_name=lambda x: 'GAMMA' if x=='G' else x
        strings='begin Kpoint_Path\n'
        for kl1,kl2 in zip(k_list,k_list[1:]):
            strings+='%-5s'%p_name(kl1[0])+' %5.2f %5.2f %5.2f  '%tuple(kl1[1])
            strings+='%-5s'%p_name(kl2[0])+' %5.2f %5.2f %5.2f\n'%tuple(kl2[1])
        strings+='end Kpoint_Path\n\n'
        return strings
    def get_projection():
        strings='begin projections\n'
        for prj in projection:
            strings+='%s:%s\n'%prj
        strings+='end projections\n\n'
        return strings
    def get_grid(k_num):
        if isinstance(k_num,int):
            klist=tuple([k_num]*3)
        elif isinstance(k_num,list):
            if len(k_num)==1:
                klist=tuple(k_num*3)
            elif len(k_num)==2:
                klist=(k_num[0],k_num[0],k_num[1])
            elif len(k_num)==3:
                klist=tuple(k_num)
            else:
                print('you should set length of k_num < 4')
                exit()
        else:
            print('k_num type is only int or list')
            exit()
        strings='mp_grid = %5d %5d %5d\n\n'%klist
        return strings
    ef=get_ef(prefix,'nscf')

    num_val=[['num_bands',nband],['num_wann',nwann],['num_iter',300]]
    if sw_so:
        num_val=num_val+[['spinors','.True.']]
    num_strings=win_strings(num_val,'num')

    dis_val=[['dis_win_max',dis_win[1]+ef],['dis_win_min',dis_win[0]+ef],
             ['dis_froz_max',frz_win[1]+ef],['dis_froz_min',frz_win[0]+ef],['dis_num_iter',300]]
    dis_strings=win_strings(dis_val,'dis')

    plot_val=[['fermi_surface_plot',TorF(sw_fs_plot)],['fermi_energy',ef],
              ['fermi_surface_num_points',fermi_mesh],['bands_plot','.True.'],['write_hr','.True.'],['write_u_matrices',TorF(uwrite)]]
    plot_strings=win_strings(plot_val,'plot')

    if sw_post_wan:
        bmu_min=ef-mu_range[0]
        bmu_max=ef+mu_range[1]
        bt_min=temp_range[0]
        bt_max=temp_range[1]
        post_val=[['boltzwann','.True.'],['boltz_kmesh',boltz_kmesh],['boltz_relax_time',btau],
                  ['boltz_mu_min',bmu_min],['boltz_mu_max',bmu_max],['boltz_mu_step',bmu_step],
                  ['boltz_temp_min',bt_min],['boltz_temp_max',bt_max],['boltz_temp_step',bt_step]]
        post_strings=win_strings(post_val,'post')
    else:
        post_strings=''
    k_path_strings=get_k_point_path()
    unit_cell_strings='begin unit_cell_cart\nbohr\n'+cell_parameter_stream(axis,deg)+'end unit_cell_cart\n\n'
    atom_strings='begin atoms_frac\n'+atom_position(atom,atomic_position)+'end atoms_frac\n\n'
    prj_strings=get_projection()
    k_grid_strings=get_grid(k_mesh_wannier)
    k_point_strings='begin kpoints\n'+k_cube_stream(k_mesh_wannier,F,T)+'end kpoints\n'

    fname='%s.win'%prefix
    fstream=(num_strings+dis_strings+plot_strings+post_strings+k_path_strings+unit_cell_strings
             +atom_strings+prj_strings+k_grid_strings+k_point_strings)
    write_file(fname,fstream)

def make_ph_in():
    fname='%s.ph'%prefix
    fstream='%s\n'%prefix
    var_inputph=['tr2_ph','alpha_mix','niter_ph','prefix','fildyn','trans',
                 'ldisp','lqdir','recover','outdir']
    if ithreads_num==0:
        var_inputph+=['fildvscf']
    if sw_ep:
        if occupations=='tetrahedra_opt':
            ep_setting="''" #initial setting in calc. of e-p int. using opt_tetrahedrn
            if sw_qshift:
                var_inputph+=['lshift_q']
        else:
            ep_setting="'interpolated'"
            var_inputph+=['electron_phonon']
    else:
        ep_setting=''
    if start_q!=0:
        var_inputph+=['start_q']
    if last_q!=0 and start_q<=last_q:
        var_inputph+=['last_q']
    if start_q==1 and last_q==1 and mpol:
        var_inputph+=['epsil']
    var_inputph+=['nq1','nq2','nq3']
    val_inputph={'tr2_ph':w_conv(ph_conv),'prefix':"'%s'"%prefix,'fildyn':"'%s'"%fildyn,'fildvscf':fildvscf,
                 'outdir':"'%s'"%outdir,'trans':'.True.','ldisp':'.True.','lqdir':'.True.','lshift_q':'.True.',
                 'recover':recover,'electron_phonon':ep_setting,'alpha_mix':amix,'niter_ph':maxiter_ph,
                 'epsil':'.True.','start_q':start_q,'last_q':last_q,
                 'nq1':q_mesh_dyn[0],'nq2':q_mesh_dyn[1],'nq3':q_mesh_dyn[2]}
    fs_inputph=make_fstring_obj('inputph',var_inputph,val_inputph,'ph')
    fstream+=fs_inputph
    write_file(fname,fstream)

def make_q2r():
    fname='%s.q2r'%(prefix)
    fstream=''
    inputq2r=['fildyn','flfrc']
    if sw_ep:
        inputq2r+=['la2F']
    fs_input=make_fstring_obj('input',inputq2r,{'fildyn':"'%s'"%(fildyn+dxml),'flfrc':"'%s'"%flfrc,'la2F':'.True.'},'q2r')
    fstream+=fs_input
    write_file(fname,fstream)

def make_matdyn(phband):
    asr="'crystal'"
    input_list=['flfrc','asr','dos']
    if phband:
        fext='freq'
        dos=False
        input_val={'flfrc':"'%s'"%(flfrc+dxml),'asr':asr,'dos':'.False.','la2F':'.True.'}
    else:
        fext='matdyn'
        dos=True
        nk1=q_mesh_dos
        nk2=q_mesh_dos
        nk3=q_mesh_dos
        input_list=input_list+['nk1','nk2','nk3']
        input_val={'flfrc':"'%s'"%(flfrc+dxml),'asr':asr,'dos':'.True.','la2F':'.True.','nk1':nk1,'nk2':nk2,'nk3':nk3}
    if sw_ep:
        input_list=input_list+['la2F']
    fname='%s.%s'%(prefix,fext)
    fstream=''
    fs_input=make_fstring_obj('input',input_list,input_val,'matdyn')
    fstream+=fs_input
    if phband:
        mat=get_cr_mat(num_brav)
        avec=mat*axis
        alat=(np.sqrt(sum(avec[0]**2)) if ibrav==0 else axis[0])
        q_list=[[kl[0],list(np.linalg.inv(mat*axis/alat).dot(kl[1]))] for kl in k_list]
        fs_qlist=k_line_stream(q_mesh_bands,q_list)
        fstream+=fs_qlist
    write_file(fname,fstream)

def make_plotband_in(mode,win):
    check_type(mode,bool)
    if mode:
        fname='eband.plotband'
        name=prefix+'_bands'
        name_out=prefix
        fext='dat'
        ef=get_ef(prefix,'scf')
        win_min=win[0]+ef
        win_max=win[1]+ef
    else:
        fname='pband.plotband'
        name='matdyn'
        name_out='matdyn'
        fext='freq'
        ef=0.
        win_min=0
        win_max=win[1]
    erange=win_max-win_min
    if 1000<erange:
        egrid=400
    elif 500<erange<=1000:
        egrid=200
    elif 100<erange<=500:
        egrid=100
    elif 50<erange<=100:
        egrid=20
    elif 10<erange<=50:
        egrid=10
    elif 5<erange<=10:
        egrid=2
    elif 1<erange<=5:
        egrid=1
    else:
        egrid=0.2
    format=(name,fext,win_min,win_max,name_out,name_out,ef,egrid,ef)
    fstream='%s.%s\n%7.3f %7.3f\n%s.xmgr\n%s.ps\n%9.5f\n%6.2f %6.2f\n'%format
    write_file(fname,fstream)

def make_epw():
    fname='%s.epw'%prefix
    dvscfdir='../save'
    val_epw={'prefix':"'%s'"%prefix,'outdir':"'./'",'dvscf_dir':"'%s'"%dvscfdir,
             'ep_coupling':'.True.','elph':'.True.',
             'nk1':k_mesh_wannier[0],'nk2':k_mesh_wannier[1],'nk3':k_mesh_wannier[2],
             'nq1':q_mesh_dyn[0],'nq2':q_mesh_dyn[1],'nq3':q_mesh_dyn[2],
             'nkf1':epw_k_mesh[0],'nkf2':epw_k_mesh[1],'nkf3':epw_k_mesh[2],
             'nqf1':epw_q_mesh[0],'nqf2':epw_q_mesh[1],'nqf3':epw_q_mesh[2]}
    var_epw=['prefix']
    for i,atm in enumerate(atom):
        tmp='amass(%d)'%(i+1)
        var_epw+=[tmp]
        val_epw.update({tmp:mass[atm]})
    var_epw+=(['outdir','ep_coupling','elph','dvscf_dir']+['nk%d'%(i+1) for i in range(3)]+['nq%d'%(i+1) for i in range(3)]
             +['nkf%d'%(i+1) for i in range(3)]+['nqf%d'%(i+1) for i in range(3)])
    fstream=make_fstring_obj('inputepw',var_epw,val_epw,'epw')
    write_file(fname,fstream)

def make_cppp():
    fname='%s.cppp'%prefix

#---------------------------- main ---------------------------------
def main(prefix):
    date()
    mpiexe=mpiopt+'mpirun -np %d '%mpi_num if sw_mpi else ''
    npool='-nk %d '%kthreads_num if kthreads_num!=0 else ''
    nimage='-ni %d '%ithreads_num if ithreads_num!=0 else ''
    #wan_exe='wan.x'
    wan_exe='wannier_ham.x' #for QE6.5 or greater
    def os_io(prefix,exe):
        name='%s.%s'%(prefix,exe)
        return '<%s>%s.out'%tuple([name]*2)
    if sw_opt: #optimaization
        make_pw_cp_in('vc-relax',False) #make pw.x's input file for scf
        if sw_run:
            ck_conv=True
            os_and_print(mpiexe+'pw.x '+npool+os_io(prefix,'scf'))
            date()
            for i in range(10):
                fname='%s.%s.out'%(name,ext)
                if os.path.exists(fname):
                    for f in open(fname,'r'):
                        if f.find('convergence NOT achieved after 100 iterations: stopping')!=-1:
                            restart="'restart'"
                            ck_conv=False
                else:
                    print('optimization error!')
                if(not ck_conv):
                    make_pw_cp_in('vc-relax',False,restart) #make pw.x's input file for scf
                    os_and_print(mpiexe+'pw.x '+npool+os_io(prefix,'scf'))
                    date()
    elif sw_md:
        make_pw_cp_in('vc-md' if sw_vc else 'md',False)
        if sw_run:
            os_and_print(mpiexe+'pw.x '+npool+os_io(prefix,'md'))
            date()
    else:
        if sw_scf: #calculate scf cycle
            make_pw_cp_in('scf',False) #make pw.x's input file for scf
            if sw_vdW and vdW_corr=='vdW-DF' and (not os.path.isfile('vdW_kernel_table')):
                os_and_print(mpiexe+'generate_vdW_kernel_table.x ')
            if sw_run:
                os_and_print(mpiexe+'pw.x '+npool+os_io(prefix,'scf'))
                date()
        if sw_dos:
            make_pw_cp_in('nscf',False) #make pw.x's input file for nscf
            if sw_prj:
                make_prjwfc_in()
            else:
                make_dos_in()
            if sw_run:
                os_and_print(mpiexe+'pw.x '+npool+os_io(prefix,'nscf'))
                if sw_prj:
                    os_and_print(mpiexe+'projwfc.x '+npool+os_io(prefix,'prjwfc'))
                else:
                    os_and_print('dos.x '+os_io(prefix,'dos_in'))
                date()
        if sw_bands: #calculate band structure
            make_pw_cp_in('bands',True) #make pw.x's input file for bands
            if sw_run:
                os_and_print(mpiexe+'pw.x '+npool+os_io(prefix,'bands'))
                date()
            make_bands_in() #make input file for bands.x
            if sw_run:
                os_and_print('bands.x '+os_io(prefix,'bands_in'))
            make_plotband_in(True,eband_win) #make input file for plotband.x
            if sw_run:
                os_and_print('plotband.x '+os_io('eband','plotband'))
        if sw_wan: #wannierize
            if sw_wan_init_nscf: #calculate energy of all k-points for wannier90
                make_pw_cp_in('nscf',True) #make pw.x's input file for nscf
                if sw_run:
                    os_and_print(mpiexe+'pw.x '+npool+os_io(prefix,'nscf'))
                    date()
            if sw_wan_init: #execute wannier90.x with initial settings
                make_win() #make input file for wannier90 (and pw2wan)
                if sw_run:
                    os_and_print('%s -pp %s'%(wan_exe,prefix))
                make_pw2wan_in()
                if sw_run:
                    os_and_print(mpiexe+'pw2wannier90.x '+os_io(prefix,'pw2wan'))
                    os_and_print('%s %s'%(wan_exe,prefix))
                    date()
            else: #execute only wannier90
                make_win()
                os_and_print('%s %s'%(wan_exe,prefix))
        if sw_ph: #calculate phonon
            make_ph_in() #make input file for ph.x
            if sw_run:
                print(mpiexe+'ph.x '+nimage+npool+os_io(prefix,'ph'))
                os_and_print(mpiexe+'ph.x '+nimage+npool+os_io(prefix,'ph'))
                date()
        if sw_dyn:
            make_q2r() #make input file for q2r.x
            if sw_run:
                if sw_so:
                    os_and_print('cp %s.dyn0 %s.dyn0.xml'%(prefix,prefix))
                os_and_print(mpiexe+'q2r.x '+npool+os_io(prefix,'q2r'))
            make_matdyn(False) #make input file for matdyn.x (for dos and alpha^2F calculation)
            if sw_run:
                os_and_print(mpiexe+'matdyn.x '+npool+os_io(prefix,'matdyn'))
                date()
                if sw_gen_a2f: #create a2F.dos.dat files for plotting
                    for i in range(qnum):
                        fname='a2F.dos%s'%(i+1)
                        tmp1=[f for f in open(fname,'r')]
                        tmp=tmp1[5:-1]
                        tmp1=[t1.strip()+t2 for t1,t2 in zip(tmp,tmp[1:])]
                        tmp=tmp1[::2]
                        f=open(fname+'.dat','w')
                        for i in tmp:
                            f.write(i)
                        f.close()
            make_matdyn(True) #make input file for matdyn.x (for phonon dispersion calculation)
            if sw_run:
                os_and_print(mpiexe+'matdyn.x '+npool+os_io(prefix,'freq'))
            make_plotband_in(False,pband_win) #make input file for plotband.x
            if sw_run:
                os_and_print('plotband.x '+os_io('pband','plotband'))
        if sw_save_dir:
            import shutil
            f=open('%s.dyn0'%prefix)
            f.readline()
            nq=int(f.readline())
            f.close()
            if not os.path.isdir('save'):
                os.makedirs('save',exist_ok=True)
            if not os.path.isdir('save/'+prefix+'.phsave'):
                shutil.copytree('_ph0/'+prefix+'.phsave','save/'+prefix+'.phsave')
            for i in range(nq):
                qn=str(i+1)
                dvscf_dir='' if i==0 else prefix+'.q_'+qn+'/'
                shutil.copy(prefix+'.dyn'+qn+dxml,'save/'+prefix+'.dyn_q'+qn)
                shutil.copy('_ph0/'+dvscf_dir+prefix+'.dvscf1','save/'+prefix+'.dvscf_q'+qn)
                shutil.copy('_ph0/'+dvscf_dir+prefix+'.dvscf_paw1','save/'+prefix+'.dvscf_paw_q'+qn)
        if sw_epw:
            make_epw()
            if sw_run:
                os_and_print(mpiexe+'epw.x '+npool+os_io(prefix,'epw'))
                date()

if __name__=="__main__":
    try: #detect k_list
        k_list
    except NameError:
        ibrav,num_brav=get_bravs(space,ibrav)
        k_list=generate_klist(num_brav)
    if sw_sc:
        axis=axis*np.array(sc_size)
        atomic_position=gen_SC_positions(sc_size,atomic_position)
        try:
            imp_atom
            if len(imp_atom)!=0:
                atom=imp_atom+atom
                try:
                    imp_pot_type
                    pot_type=imp_pot_type+pot_type
                except NameError:
                    print('conflict atomic num')
                    exit()
        except NameError:
            print('no impurity',flush=True)
        if sc_size[0]!=sc_size[1]:
            if num_brav in {1,2,3,4,5,6,7}:
                space=1
        elif sc_size[0]!=sc_size[2]:
            if num_brav in {1,2,3,5}:
                space=1
    if from_poscar:
        axis,cry_ax,atom,atomic_position=read_poscar()
        ibrav=0
    ibrav,num_brav=get_bravs(space,ibrav)
    try: #detect type_xc
        type_xc
    except NameError:
        type_xc='pbe'
    if sw_so or sw_apw:
        txc='rel-'+type_xc if sw_so else type_xc
        pseude_dir='/home/suzu/pslibrary/%s/PSEUDOPOTENTIALS/'%txc
        #pseude_dir='/home/usr2/h70252j/UPF/%s/'%txc
    else:
        txc=type_xc
        pseude_dir='/home/suzu/pslibrary/%s/PSEUDOPOTENTIALS/'%txc
        #pseude_dir='/home/usr2/h70252j/UPF/'
    UPF=['%s.%s-'%(at[:2],txc)+pp+'-'+pot_kind+'.UPF' for pp, at in zip(pot_type,atom)]
    main(prefix)

#==============================================================================#
# This is auto calculation script for Quantum Espresso                         #
#                                                                              #
# Copyright (c) 2018-2021  K. Suzuki                                           #
#==============================================================================#
#                          The MIT License (MIT)                               #
#==============================================================================#
#Permission is hereby granted, free of charge, to any person obtaining a copy  #
#of this software and associated documentation files (the "Software"), to deal #
#in the Software without restriction, including without limitation the rights  #
#to use, copy, modify, merge, publish, distribute, sublicense, and/or sell     #
#copies of the Software, and to permit persons to whom the Software is         #
#furnished to do so, subject to the following conditions:                      #
#                                                                              #
#The above copyright notice and this permission notice shall be included in all#
#copies or substantial portions of the Software.                               #
#                                                                              #
#THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR    #
#IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,      #
#FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE   #
#AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER        #
#LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, #
#OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE #
#SOFTWARE.                                                                     #
#==============================================================================#
