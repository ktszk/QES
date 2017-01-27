#!/usr/bin/env python
# -*- coding:utf-8
#$ -cwd
#$ -V -S /usr/local/Python-2.7.12/bin/python
#$ -N NbSe2
#$ -e out.e
#$ -o out.o
#$ -pe fillup 64
#$ -q lantana.q
#BSUB -q "queue_name"
#BSUB -n 16
#BSUB -J "Job name"
#BSUB -o std_output_file
#BSUB -m "node_name"
(T,F)=(True,False)                #alias for bool numbers 
(mpi_num, kthreads_num) = (64, 8) #number of threads for MPI/openMP

#optional parameters
aa=3.4446   #lattice parameter a
ab=aa       #lattice parameter b
ac=12.5444  #lattice parameter c
za=0.1172   #non-toribial atomic position

#======================Crystal structure===========================
prefix='NbSe2'                          #material name (or job name)
space=194                               #space group
axis=[aa,ab,ac]                         #lattice parameters a,b,c
deg=[90,90,120]                         #lattice parameters alpha,beta,gamma
atom=['Nb','Se']                        #elements name
atomic_position=[[[0. ,0. ,0.25], [0. ,0. ,0.75]],
                 [[1./3 ,2./3 ,za], [2./3 ,1./3 ,1-za], 
                  [1./3 ,2./3 ,0.5-za], [2./3 ,1./3 ,0.5+za]]]
#------------------------------------------------------------------
ibrav=0                                 #brave lattice type
type_xc='pbe'                           #type of exchange correlation functional
pot_type=['%s.%s-nsp-van',
          '%s.%s-n-kjpaw_psl.0.2']      #psede potential name
#====================directorys settings===========================
pseude_dir='/home/Apps/upf_files/'      #path of psede potential directory
outdir='./'                             #path of output directory
#===============switch & number of parallel threads================
sw_scf = T                              #generate input file for scf calculation
sw_bands = T                            #generate input file for band calculation
sw_ph = T                               #generate input file for phonon calculation
sw_wan = T                              #switch wannierization
sw_wan_init = F                         #generate input file for wannier90
sw_wan_init_nscf = F                    #calc nscf cycle for wannier90
sw_restart = F                          #switch restart tag or not

sw_run = F                              #switch of execute DFT calculation or not
sw_mpi = T                              #switch of MPI calculation
sw_bsub = F                             #switch bsub
#=====================pw_parameters================================
k_mesh_scf=[16, 16, 8]                  #k mesh for DFT calc
k_mesh_bands=10                         #k mesh for bands calc
k_mesh_wannier=8                        #k mesh for wannierize
(ecut, ec_rho)=(60.0, 800)              #cut off energy of pw and density
(e_conv, f_conv)=(1.0e-5, 1.0e-4)       #threshold of total energy's convergence and force's one 
(scf_conv,nscf_conv)=(1.0e-12, 1.0e-10) #threshold of convergence on scf,nscf cycles
nband=50                                #number of bands
nstep=500                               #max number of scf cycle's step
deg=0.025                               #dispersion of k-mesh
eband_win=[-3, 3]                       #energy range of .ps file
#======================ph_parameters===============================
q_mesh_dyn=[8, 8, 2]                    #q mesh for phonon DFPT calc
q_mesh_bands=20                         #q mesh for phonon band calc
q_mesh_dos=8                            #q mesh for phonon dos calc
ph_conv=1.0e-14                         #threshold of energy's convergence for phonon
pband_win=[0, 350]                      #energy range of .ps file
#===================Wannier_parameters=============================
nwann=4                                 #number of wannier basis
dis_win=[-2.80, 1.40]                   #max(min)_window, range of sub space energy
frz_win=[-0.60, 1.00]                   #froz_window
projection=[('Ru', 'dxz,dyz,dxy')]      #projections, initial funcution of wannier
sw_fs_plot= F                           #plot Fermi surface
fermi_mesh = 100                        #mesh of k-points in bxsf file
unk=F                                   #Bloch(Wannier)_func
uwrite=T                                #output unitary matrix for bloch to wannier
#=============================modules==============================
import numpy as np
import os, datetime
#=====================global_lambda_expression=====================
TorF=lambda x:'.True.' if x else '.False.'
w_conv=lambda a:'1.0E%d'%int(np.log10(a))
#========== list of k-point for band calculation (manual) =========
#k_list=[['G',[0.,0.,0.]],['M',[0.5,0.,0.]],['K',[1./3,1./3,0.]],
#        ['G',[0.,0.,0.]],['Z',[0.,0.,0.5]]] #Phexa
#====================== physical parameters =======================
bohr=round(0.52917721067, 6)           #Bohr Radius
ibohr=1.0/bohr                         #inverse of Bohr Radius
#========================= atomic mass ============================
mass={'H':1.00794,                                                                            'He':4.002602,
      'Li':6.941,'Be':9.012182,'B':10.811,'C':12.0107, 'N':14.0067,'O':15.9994,'F':18.9984032,'Ne':20.1797,
      'Na':22.98977,'Mg':24.305,  'Al':26.981538,'Si':28.0855,'P':30.973761,'S':32.065,'Cl':35.453,'Ar':39.948,
      'K':39.0983,  'Ca':40.078,  'Sc':44.95591,
      'Ti':47.867,'V':50.9415,'Cr':51.9961,'Mn':54.938049,'Fe':55.845,'Co':58.9332,'Ni':58.6934,'Cu':63.546,
      'Zn':65.409,'Ga':69.723,'Ge':72.64,'As':74.9216,'Se':78.96,'Br':79.904,'Kr':83.798,
      'Rb':85.4678,'Sr':87.62,'Y':88.90585,
      'Zr':91.224,'Nb':92.90638,'Mo':95.94,'Tc':97.907216,'Ru':101.07,'Rh':102.9055,'Pd':106.42,'Ag':107.8682,
      'Cd':112.411,'In':114.818,'Sn':118.71,'Sb':121.76,'Te':127.6,'I':126.90447,'Xe':131.293,
      'Cs':132.90545,'Ba':137.327,
      'La':138.9055,'Ce':140.116,'Pr':140.90765,'Nd':144.24,'Pm':144.912744,'Sm':150.36,'Eu':151.964,'Gd':157.25,
      'Tb':158.92534,'Dy':162.5,'Ho':164.93032,'Er':167.259,'Tm':168.93421,'Yb':173.04,'Lu':174.967,
      'Hf':178.49,'Ta':180.9479,'W':183.84,'Re':186.207,'Os':190.23,'Ir':192.217,'Pt':195.078,'Au':196.96655,
      'Hg':200.59,'Tl':204.3833,'Pb':207.2,'Bi':208.98038,
      'Th':232.0381,'Pa':231.03588,'U':238.02891}
#==========================global_variables========================
UPF=[pp%(at,type_xc)+'.UPF' for pp, at in zip(pot_type,atom)]
axis=np.array(axis)
fildvscf="'dvscf'"
fildyn="'%s.dyn'"%prefix
flfrc="'%s.fc'"%prefix
recover=TorF(sw_restart)
restart="'from_scratch'"
#==============serch several variables and initialize if can't find
try: #detect type_xc
    type_xc
except NameError:
    type_xc='pbe'

try: #detect brav
    brav
    try:
        hexa
    except NameError:
        if isinstance(space, int):
            hexa=T if space in range(168,195) else F
        elif isinstance(space, str):
            hexa=T if '6' in space else F
except NameError:
    hexa=F
    if isinstance(space,int):
        if space in {23,24,44,45,46,71,72,73,74,79,80,82,87,88,
                     97,98,107,108,109,110,119,120,121,122,139,
                     140,141,142,197,199,204,206,211,214,217,220,229,230}:
            brav='I'
        elif space in {22,42,43,69,70,196,202,203,209,210,216,219,225,226,227,228}:
            brav='F'
        elif space in {146,148,155,160,161,166,167}:
            brav='R'
        elif space in {38,39,40,41}:
            brav='A'
        elif space in {5,8,9,12,15,20,21,35,36,37,63,64,65,66,67,68}:
            brav='C'
        else:
            brav='P'
            if space in range(168,195): #168>194 is Hexagonal
                hexa=T
    elif isinstance(space,str):
        if 'I' in space:
            brav='I'
        elif 'F' in space:
            brav='F'
        elif ('R' in space):
            brav='R'
        elif ('A' in space):
            brav='A'
        elif ('C' in space):
            brav='C'
        else:
            brav='P'
            if '6' in space:
                hexa=T

try: #detect k_list
    k_list
except NameError: #common k-points list
    if brav=='I':
        if axis[0]==axis[1]:
            if axis[0]==axis[2]: #cube
                k_list=[['G',[0.,0.,0.]],['H',[0.5,0.5,0.5]],
                        ['P',[0.25,0.75,-0.25]],['N',[0.,0.5,0.]],
                        ['G',[0.,0.,0.]]]
            else: #tetra
                k_list=[['G',[0.,0.,0.]],['X',[0.0,0.5,-0.5]],
                        ['P',[0.25,0.75,-0.25]],['N',[0.,0.5,0.]],
                        ['G',[0.,0.,0.]],['Z',[0.5,0.5,0.5]]]
    elif brav=='F':
        if axis[0]==axis[2]: #cube
            k_list=[['G',[0.,0.,0.]],['X',[0.5,0.,0.5]],['G',[1.0,0.,0.]]
                    ['L',[0.5,0.5,0.5]],['W',[0.5,0.25,0.75]],
                    ['G',[0.,0.,0.]]]
        else: #ortho
            pass
    else:
        if hexa: #Phexa
            k_list=[['G',[0.,0.,0.]],['K',[2./3,-1./3,0.]],
                    ['M',[0.5,0.,0.]],['G',[0.,0.,0.]],['Z',[0.,0.,0.5]]]
            #k_list=[['G',[0.,0.,0.]],['K',[2./3,0.,0.]],
            #        ['M',[0.5,-0.5/np.tan(2.*np.pi/3.),0.]],
            #        ['G',[0.,0.,0.]],['Z',[0.,0.,0.5]]]
        else:
            if axis[0]==axis[1]:
                if axis[0]==axis[2]: #cube
                    k_list=[['R',[0.5,0.5,0.]],['G',[0.,0.,0.]]
                            ,['X',[0.5,0.0,0.]],['M',[0.5,0.5,0.]],
                            ['G',[0.,0.,0.]]]
                else: #tetra
                    k_list=[['G',[0.,0.,0.]],['X',[0.5,0.0,0.]],
                            ['M',[0.5,0.5,0.]],['G',[0.,0.,0.]],
                            ['Z',[0.,0.,0.5]]]
if sw_fs_plot:
    try: #detect fermi_mesh
        fermi_mesh
    except NameError:
        fermi_mesh=100
#==========================functions===============================
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
                     else '%a %b %d %X JST %Y'))

def os_and_print(command):
    print(command)
    os.system(command)

def get_ef(name,ext):
    fname='%s.%s.out'%(name,ext)
    if os.path.exists(fname):
        for f in open(fname,'r'):
            if f.find('Fermi')!=-1:
                item=f.split('is')
                it1=item[1].split('ev')
                return float(it1[0])
        else:
            print('input error from %s.%s.out\n'%(name,ext)
                  +'return ef = 0\n')
            return 0.
    else:
        print('can not find %s\n return ef = 0\n'%fname)
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
        fstring=fstring+form%vc+str(val_dic[vc])+'\n'
    fstring=fstring+'/\n'
    return fstring

def atom_position(atom,atomic_position):
    mat=get_cr_mat(brav,False)
    mat=np.linalg.inv(mat).T
    aposition=[[list(mat.dot(np.array(ap))) for ap in app] for app in atomic_position]
    atom_string=''
    for i,at in enumerate(atom):
        for ap in aposition[i]:
            atom_string=atom_string+'%2s  %12.9f %12.9f %12.9f\n'%tuple([at]+ap)

    return atom_string

def get_cr_mat(brav,hexa):
    if brav=='P':
        if hexa:
            mat=np.array([[ 1. ,0.             ,0.],
                          [-0.5,np.sqrt(3.)*0.5,0.],
                          [ 0. ,0.             ,1.]])
        else:
            mat=np.identity(3)
    elif brav=='I':
        mat=np.array([[ 0.5,-0.5,0.5],
                      [ 0.5, 0.5,0.5],
                      [-0.5,-0.5,0.5]])
    elif brav=='F':
        mat=np.array([[-0.5,0. ,0.5],
                      [ 0. ,0.5,0.5],
                      [-0.5,0.5,0.]])
    elif brav=='C':
        mat=np.array([[ 0.5,0.5,0.],
                      [-0.5,0.5,0.],
                      [ 0. ,0. ,1.]])
    elif brav=='R':
        mat=np.array([[ 1.0,0.,0.],
                      [np.cos(np.pi*deg[0]/180),np.cos(np.pi*deg[0]/180),0.],
                      [np.cos(np.pi*deg[0]/180),0. ,1.]])
    elif brav=='A':
        mat=np.identity(3)
    return mat

def cell_parameter_stream(axis,deg):
    mat_ax=np.identity(3)*axis*ibohr
    mat=get_cr_mat(brav,hexa)
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
            k_string=k_string+'%9.6f %9.6f %9.6f'%(ks[0]+kp[0]*i,ks[1]+kp[1]*i,ks[2]+kp[2]*i)+wt+'\n'
    k_string=k_string+'%9.6f %9.6f %9.6f'%tuple(k_list[-1][1])+wt+'\n'
    return k_string

def k_cube_stream(k_num,w_sw,sw_wan):
    wfunc=lambda x,y,z:'\n' if x else ' %f\n'%(1./z if y else 1.)
    if not isinstance(w_sw,bool):
        w_sw=False
    if isinstance(k_num,int):
        dk=1./k_num
        wst=wfunc(sw_wan,w_sw,(k_num**3))
        k_string='' if sw_wan else '%d\n'%k_num**3
        for i in range(k_num):
            for j in range(k_num):
                for k in range(k_num):
                    k_string=k_string+'%f %f %f'%(dk*i,dk*j,dk*k)+wst
    elif isinstance(k_num,list):
        if  len(k_num==2):
            dk=[1./kp for kp in k_num]
            dk=[dk[0]]+dk
            kn=[k_num[0]]+k_num
        elif len(k_num==3):
            dk=[1./kp for kp in k_num]
            kn=k_num
        else:
            print('k dimension < 3')
            exit()
        if w_sw:
            w0=1.
            for wi in range(kn):
                w0=w0*wi
        wst=wfunc(sw_wan,w_sw,w0)
        k_string='' if sw_wan else '%d\n'%w0
        for i in range(kn[0]):
            for j in range(kn[1]):
                for k in range(kn[2]):
                    k_string=k_string+'%f %f %f'%(dk[0]*i,dk[1]*j,dk[2]*k)+wst
    else:
        print('Please input list or int into k_point ')
        exit()

    return k_string
#---------------------input file generators-------------------------------
def make_pw_in(calc):
    def atomic_parameters_stream(atom,atomic_position,UPF):
        atom_string='\nATOMIC_SPECIES\n'
        for at,up in zip(atom,UPF):
            atom_string=atom_string+' %-2s %11.7f  %s\n'%(at,mass[at],up)
        atom_string=atom_string+'\nATOMIC_POSITIONS crystal\n'
        atom_string=atom_string+atom_position(atom,atomic_position)
        atom_string=atom_string+'\n'
        return atom_string

    (fext,convthr) =(('nscf' if calc=='nscf' else 'bands',nscf_conv) 
                     if calc in ['nscf','bands'] else ('scf',scf_conv))
    fname='%s.%s'%(prefix,fext)
    fstream=''
    var_control=['title','calculation','restart_mode','outdir','pseudo_dir',
                 'prefix','etot_conv_thr','forc_conv_thr','nstep','tstress','tprnfor','wf_collect']
    val_control={'title':"'%s'"%prefix,'calculation':"'%s'"%calc,'restart_mode':restart,'outdir':"'%s'"%outdir,
                 'pseudo_dir':"'%s'"%pseude_dir,'prefix':"'%s'"%prefix,'etot_conv_thr':w_conv(e_conv),
                 'forc_conv_thr':w_conv(f_conv),'nstep':nstep,'tstress':'.True.','tprnfor':'.True.',
                 'wf_collect':'.True.'}
    fs_control=make_fstring_obj('control',var_control,val_control,'pw')
    fstream=fstream+fs_control

    var_system=['ibrav','nat','ntyp','occupations','smearing','degauss',
                'la2f','nbnd','ecutwfc','ecutrho']
    val_system={'ibrav':ibrav,'nat':sum(len(a) for a in atomic_position),
                'ntyp':len(atom),'occupations':"'smearing'",'smearing':"'marzari-vanderbilt'",
                'degauss':0.025,'la2f':'.True.','nbnd':nband,'ecutwfc':ecut,'ecutrho':ec_rho}
    fs_system=make_fstring_obj('system',var_system,val_system,'pw')
    fstream=fstream+fs_system

    var_electrons=['diagonalization','conv_thr']
    val_electrons={'diagonalization':"'david'",'conv_thr':w_conv(convthr)}
    fs_electrons=make_fstring_obj('electrons',var_electrons,val_electrons,'pw')
    fstream=fstream+fs_electrons

    if calc in {'relax','md','vc-relax','vc-md'}:
        var_ions=[]
        val_ions={}
        fs_ions=make_fstring_obj('ions',var_ions,val_ions,'pw')
        fstream=fstream+fs_ions

    if calc in {'vc-relax','vc-md'}:
        var_cell=[]
        val_cell={}
        fs_cell=make_fstring_obj('cell',var_cell,val_cell,'pw')
        fstream=fstream+fs_cell

    fs_atomparam=atomic_parameters_stream(atom,atomic_position,UPF)
    fstream=fstream+fs_atomparam
    if ibrav==0:
        fs_cellparam='CELL_PARAMETERS\n'+cell_parameter_stream(axis,deg)+'\n'
        fstream=fstream+fs_cellparam

    fstream=fstream+'K_POINTS (%s)\n'%('crystal' if calc in {'nscf', 'bands'} else 'automatic')
    if calc in {'nscf', 'bands'}:
        if calc=='nscf':
            fs_kl_param=k_cube_stream(k_mesh_wannier,T,F)
        else:
            fs_kl_param=k_line_stream(k_mesh_bands,k_list)
        fstream=fstream+fs_kl_param
    else:
        if len(k_mesh_scf)==3:
            k_mesh=k_mesh_scf
        elif len(k_mesh_scf)==2:
            k_mesh=[k_mesh_scf[0]]*2+[k_mesh_scf[1]]
        else:
            k_mesh=[k_mesh_scf]*3
        fstream=fstream+'%d %d %d %d %d %d\n'%tuple(k_mesh+[0]*3)
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
    fstream=fstream+fs_bands
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
            strings=strings+form%vs[0]+format_val(vs[1])+'\n'
        strings=strings+'\n'
        return strings
    def get_k_point_path():
        p_name=lambda x: 'GAMMA' if x=='G' else x
        strings='begin Kpoint_Path\n'
        for kl1,kl2 in zip(k_list,k_list[1:]):
            strings=strings+'%-5s'%p_name(kl1[0])+' %5.2f %5.2f %5.2f  '%tuple(kl1[1])
            strings=strings+'%-5s'%p_name(kl2[0])+' %5.2f %5.2f %5.2f\n'%tuple(kl2[1])
        strings=strings+'end Kpoint_Path\n\n'
        return strings
    def get_projection():
        strings='begin projections\n'
        for prj in projection:
            strings=strings+'%s:%s\n'%prj
        strings=strings+'end projections\n\n'
        return strings
    def get_grid(k_num):
        if isinstance(k_num,int):
            klist=tuple([k_num]*3)
        elif isinstance(k_num,int):
            if len(k_num)==2:
                klist=(k_num[0],k_num[0],k_num[1])
            else:
                klist=tuple(k_num)
        strings='mp_grid = %5d %5d %5d\n\n'%klist
        return strings
    ef=get_ef(prefix,'nscf')

    num_val=[['num_bands',nband],['num_wann',nwann],['num_iter',300]]
    num_strings=win_strings(num_val,'num')

    dis_val=[['dis_win_max',dis_win[1]+ef],['dis_win_min',dis_win[0]+ef],
             ['dis_froz_max',frz_win[1]+ef],['dis_froz_min',frz_win[0]+ef],['dis_num_iter',300]]
    dis_strings=win_strings(dis_val,'dis')

    plot_val=[['fermi_surface_plot',TorF(sw_fs_plot)],['fermi_energy',ef],
              ['fermi_surface_num_points',fermi_mesh],['bands_plot','.True.'],['write_hr','.True.'],['write_u_matrices',TorF(uwrite)]]
    plot_strings=win_strings(plot_val,'plot')

    k_path_strings=get_k_point_path()
    unit_cell_strings='begin unit_cell_cart\nbohr\n'+cell_parameter_stream(axis,deg)+'end unit_cell_cart\n\n'
    atom_strings='begin atoms_frac\n'+atom_position(atom,atomic_position)+'end atoms_frac\n\n'
    prj_strings=get_projection()
    k_grid_strings=get_grid(k_mesh_wannier)
    k_point_strings='begin kpoints\n'+k_cube_stream(k_mesh_wannier,F,T)+'end kpoints\n'

    fname='%s.win'%prefix
    fstream=(num_strings+dis_strings+plot_strings+k_path_strings+unit_cell_strings
             +atom_strings+prj_strings+k_grid_strings+k_point_strings)
    write_file(fname,fstream)

def make_ph_in():
    fname='%s.ph'%prefix
    fstream='%s\n'%prefix
    var_inputph=['tr2_ph','prefix','fildyn','trans','ldisp','recover',
                 'outdir','fildvscf','electron_phonon','nq1','nq2','nq3']
    val_inputph={'tr2_ph':w_conv(ph_conv),'prefix':"'%s'"%prefix,'fildyn':fildyn,'fildvscf':fildvscf,
                 'outdir':"'%s'"%outdir,'trans':'.True.','ldisp':'.True.','recover':recover,
                 'electron_phonon':"'interpolated'",
                 'nq1':q_mesh_dyn[0],'nq2':q_mesh_dyn[1],'nq3':q_mesh_dyn[2]}
    fs_inputph=make_fstring_obj('inputph',var_inputph,val_inputph,'ph')
    fstream=fstream+fs_inputph
    write_file(fname,fstream)

def make_q2r():
    fname='%s.q2r'%(prefix)
    fstream=''
    fs_input=make_fstring_obj('input',['fildyn','flfrc','la2F'],
                              {'fildyn':fildyn,'flfrc':flfrc,'la2F':'.True.'},'q2r')
    fstream=fstream+fs_input
    write_file(fname,fstream)

def make_matdyn(phband):
    asr="'crystal'"
    if phband:
        fext='freq'
        dos=False
        input_list=['flfrc','asr','dos']
        input_val={'flfrc':flfrc,'asr':asr,'dos':'.False.'}
    else:
        fext='matdyn'
        dos=True
        nk1=q_mesh_dos
        nk2=q_mesh_dos
        nk3=q_mesh_dos
        input_list=['flfrc','asr','dos','la2F','nk1','nk2','nk3']
        input_val={'flfrc':flfrc,'asr':asr,'dos':'.True.','la2F':'.True.','nk1':nk1,'nk2':nk2,'nk3':nk3}
    fname='%s.%s'%(prefix,fext)
    fstream=''
    fs_input=make_fstring_obj('input',input_list,input_val,'matdyn')
    fstream=fstream+fs_input
    if phband:
        mat=get_cr_mat(brav,hexa)
        q_list=[[kl[0],list(np.linalg.inv(mat*(axis/axis[0])).dot(kl[1]))] for kl in k_list]
        fs_klist=k_line_stream(q_mesh_bands,q_list)
        fstream=fstream+fs_klist
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
    format=(name,fext,win_min,win_max,name_out,name_out,ef,nband,ef)
    fstream='%s.%s\n%d %d\n%s.xmgr\n%s.ps\n%9.5f\n%6.2f %6.2f\n'%format
    write_file(fname,fstream)
#---------------------------- main ---------------------------------
def main(prefix):
    date()
    mpiopt='$LSF_BINDIR/openmpi-' if sw_bsub else ''
    mpiexe=mpiopt+'mpirun -np %d '%mpi_num if sw_mpi else ''
    npool='-nk %d '%kthreads_num if kthreads_num!=0 else ''
    wan_exe='wan.x'
    def os_io(prefix,exe):
        name='%s.%s'%(prefix,exe)
        return '<%s>%s.out'%tuple([name]*2)
    if sw_scf: #calculate scf cycle
        make_pw_in('scf') #make pw.x's input file for scf
        if sw_run:
            os_and_print(mpiexe+'pw.x '+npool+os_io(prefix,'scf'))
            date()
    if sw_bands: #calculate band structure
        make_pw_in('bands') #make pw.x's input file for bands
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
            make_pw_in('nscf') #make pw.x's input file for nscf
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
            os_and_print(mpiexe+'ph.x '+npool+os_io(prefix,'ph'))
            date()
        make_q2r() #make input file for q2r.x
        if sw_run:
            os_and_print(mpiexe+'q2r.x '+npool+os_io(prefix,'q2r'))
        make_matdyn(False) #make input file for matdyn.x (for dos and alpha^2F calculation)
        if sw_run:
            os_and_print(mpiexe+'matdyn.x '+npool+os_io(prefix,'matdyn'))
            date()
        make_matdyn(True) #make input file for matdyn.x (for phonon dispersion calculation)
        if sw_run:
            os_and_print(mpiexe+'matdyn.x '+npool+os_io(prefix,'freq'))
        make_plotband_in(False,pband_win) #make input file for plotband.x
        if sw_run:
            os_and_print('plotband.x '+os_io('pband','plotband'))
    date()

if __name__=="__main__":
    main(prefix)
