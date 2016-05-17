#!/usr/bin/env python
# -*- coding:utf-8
prefix='Sr2RuO4'

space=139
aa=3.8603; ab=aa; ac=12.729
zp1,zp=(0.35316,0.1619)

axis=[aa,ab,ac]
deg=[90,90,90]
atom=['Sr','Ru','O']
atomic_position=[[[0.,0.,zp1],[0.,0.,1-zp1]],
                 [[0.,0.,0.]],
                 [[0.,0.5,0.],[0.5,0.,0.],[0.,0.,zp],[0.,0.,1-zp]]]
ibrav=0
type_xc='pbe'
pot_type=['%s.%s-nsp-van','%s.%s-n-van','%s.%s-kjpaw_psl.0.1']
psude_dir='/home/Apps/upf_files/'

sw_bands=False
sw_ph=True
sw_wan=False
sw_restart=True

#=====================pw_parameters================================
k_mesh_scf=10
k_mesh_bands=20
k_mesh_wannier=8
ecut=60.0
ec_rho=600
conv=1.0e-9
nband=50
deg=0.025
eband_win=[-3,3]
#======================ph_parameters===============================
q_mesh_dyn=[4,4,4]
q_mesh_bands=20
q_mesh_dos=8
pband_win=[0,900]
#======================pw_&_ph_common_parameter====================
k_list=[['Gamma',0.,0.,0.],['X',0.5,0.,0.],['M',0.5,0.5,0.],['Gamma',0.,0.,0.]]
#===================Wannier_parameters=============================
nwann=10                     #number of wannier
dis_win=[-3.00,4.00]         #max(min)_window
frz_win=[-0.0,0.0]           #froz_window
projection=[['Fe','d']]      #projections
sw_fs_plot= False            #plot Fermi surface
unk=False                    #Bloch(Wannier)_func
#================physical parameter================================
bohr=round(0.52917721092,6); ibohr=1.0/bohr
mass={'H':1.00794,'He':4.002602,
      'Li':6.941,'Be':9.012182,'B':10.811,'C':12.0107,'N':14.0067,'O':15.9994,'F':18.9984032,'Ne':20.1797,
      'Na':22.98977,'Mg':24.305,'Al':26.981538,'Si':28.0855,'P':30.973761,'S':32.065,'Cl':35.453,'Ar':39.948,
      'K':39.0983,'Ca':40.078,'Sc':44.95591,
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
#=====================global_variables=============================
UPF=[pp%(at,type_xc)+'.UPF' for pp,at in zip(pot_type,atom)]
outdir='./'
fildvscf="'.dvscf'"
fildyn="'%s.dyn'"%prefix
flfrc='%s.fc'%prefix
restart='.True.' if sw_restart else '.False.'

try:
    brav
except NameError:
    if isinstance(space,int):
        if space in [23,139]:
            brav='I'
        elif space in [22]:
            brav='F'
        elif space in [146]:
            brav='R'
        elif space in [38]:
            brav='A'
        elif space in [5,8,9,12,15]:
            brav='C'
        else:
            brav='P'
    if isinstance(space,str):
        if 'I' in space:
            brav='I'
        elif 'F' in space:
            brav='F'
        elif 'R' in space:
            brav='R'
        elif 'A' in space:
            brav='R'
        elif 'C' in space:
            brav='R'
        else:
            brav='P'
#======================modules=====================================
import numpy as np
import os
#=====================functions====================================
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

def get_ef(name,ext):
    fname='%s.%s.out'%(name,ext)
    if os.path.exists(fname):
        for f in open(fname,'r'):
            if f.find('Fermi')!=-1:
                item=f.split('is')
                it1=item[1].split('ev')
                return float(it1[0])
        else:
            print('input error from %s.%s.out\n'%(name,ext))
            print('return ef = 0\n')
            return 0.
    else:
        print('not find %s\n return ef = 0\n'%fname)
        return 0.

def make_fstring_obj(obj_name,var_list,val_dic,sw_form):
    if sw_form=='pw':
        form='%28s = '
    elif sw_form=='ph':
        form='%20s = '
    else:
        form='%12s = '
    fstring='&%s\n'%obj_name
    for vc in var_list:
        fstring=fstring+form%vc+str(val_dic[vc])+'\n'
    fstring=fstring+'/\n'
    return fstring

def get_cr_mat(brav):
    if brav=='P':
        if False:
            mat=np.array([[1.,0.,0.],[-0.5,np.sqrt(3.)*0.5,0.],[0.,0.,1.]])
        else:
            mat=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    elif brav=='I':
        mat=np.array([[0.5,-0.5,0.5],[0.5,0.5,0.5],[-0.5,-0.5,0.5]])
    elif brav=='F':
        mat=np.array([[-0.5,0.,0.5],[0.,0.5,0.5],[-0.5,0.5,0.]])
    elif brav=='C':
        mat=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    elif brav=='A':
        mat=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    elif brav=='R':
        mat=np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
    return mat

def cell_parameter_stream(axis,deg):
    mat_ax=np.array([[axis[0],0,0],[0,axis[1],0],[0,0,axis[2]]])
    mat=get_cr_mat(brav)
    a_vec=list(mat.dot(mat_ax))
    cell_string='cell_parameters\n'
    for aa in a_vec:
        cell_string=cell_string+'  %12.9f  %12.9f  %12.9f\n'%tuple(aa)
    cell_string=cell_string+'\n'

    return cell_string

def atomic_parameters_stream(atom,atomic_position,UPF):
    mat=get_cr_mat(brav)
    mat=np.linalg.inv(mat).T
    aposition=[[list(mat.dot(np.array(ap))) for ap in app] for app in atomic_position]
    atom_string='\natomic_species\n'
    for at,up in zip(atom,UPF):
        atom_string=atom_string+'%2s  %11.7f  %s\n'%(at,mass[at],up)
    atom_string=atom_string+'\natomic_positions (crystal)\n'
    for i,at in enumerate(atom):
        for ap in aposition[i]:
            atom_string=atom_string+'%2s  %12.9f %12.9f %12.9f\n'%tuple([at]+ap)
    atom_string=atom_string+'\n'
    return atom_string

def k_line_stream(k_num,k_list):
    k_string='%d\n'%k_num
    dk=1./k_num
    w=1.
    wt=' %f'%w if True else ''
    for kb,ka in zip(k_list,k_list[1:]):
        kp=[kaa-kbb for kaa,kbb in zip(ka[1:],kb[1:])]
        for i in range(k_num):
            k_string=k_string+'%f %f %f'%(kp[0]*i,kp[1]*i,kp[2]*i)+wt+'\n'
    return k_string

def k_cube_stream(k_num,w_sw):
    if not isinstance(w_sw,bool):
        w_sw=False
    if isinstance(k_num,int):
        dk=1./k_num
        w=1./(k_num**3) if w_sw else 1.
        k_string='%d\n'%k_point**3
        for i in range(k_num):
            for j in range(k_num):
                for k in range(k_num):
                    k_string=k_string+'%f %f %f %f\n'%(dk*i,dk*j,dk*k,w)
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
        w=1./w0 if w_sw else 1.
        k_string='%d\n'%w0
        for i in range(kn[0]):
            for j in range(kn[1]):
                for k in range(kn[2]):
                    k_string=k_string+'%f %f %f %f\n'%(dk[0]*i,dk[1]*j,dk[2]*k,w)
    else:
        print('Please input list or int into k_point ')
        exit()

    return k_string

def make_pw_in(calc):
    fext='nscf' if calc in ['nscf','bands'] else 'scf'
    fname='%s.%s'%(prefix,fext)
    fstream=''
    var_control=['title','calculation','restart_mode','outdir','psude_dir',
                 'prefix','etot_conv_thr','forc_conv_thr','nstep','tstress','tprnfor']
    val_control={'title':prefix,'calculation':calc,'restart_mode':restart,'outdir':"'%s'"%outdir,
                 'psude_dir':"'%s'"%psude_dir,'prefix':prefix,'etot_conv_thr':conv,'forc_conv_thr':conv,
                 'nstep':100,'tstress':'.True.','tprnfor':'.True.'}
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
    val_electrons={'diagonalization':"'david'",'conv_thr':1.0e-11}
    fs_electrons=make_fstring_obj('electrons',var_electrons,val_electrons,'pw')
    fstream=fstream+fs_electrons

    if calc in ['relax','md','vc-relax','vc-md']:
        var_ions=[]
        val_ions={}
        fs_ions=make_fstring_obj('ions',var_ions,val_ions,'pw')
        fstream=fstream+fs_ions

    if calc in ['vc-relax','vc-md']:
        var_cell=[]
        val_cell={}
        fs_cell=make_fstring_obj('cell',var_cell,val_cell,'pw')
        fstream=fstream+fs_cell

    fs_atomparam=atomic_parameters_stream(atom,atomic_position,UPF)
    fstream=fstream+fs_atomparam
    if ibrav==0:
        fs_cellparam=cell_parameter_stream(axis,deg)
        fstream=fstream+fs_cellparam

    fstream=fstream+'K_POINTS (%s)\n'%('crystal' if calc in ['nscf, bands'] else 'automatic')
    if calc in ['nscf, bands']:
        if sw_cube:
            fs_kl_param=k_cube_stream(k_mesh_wannier,w_sw)
        else:
            fs_kl_param=k_line_stream(k_mesh_bands,k_list)
        fstream=fstream+fs_kl_param
    else:
        fstream=fstream+'%d %d %d %d %d %d\n'%tuple([k_mesh_scf]*3+[0]*3)
    write_file(fname,fstream)

def make_bands_in():
    pass

def make_ph_in():
    fname='%s.ph'%prefix
    fstream='%s\n'%prefix
    var_inputph=['tr2_ph','prefix','fildyn','trans','ldisp','recover',
                 'outdir','fildvscf','electron_phonon','nq1','nq2','nq3']
    val_inputph={'tr2_ph':1.0e-14,'prefix':"'%s'"%prefix,'fildyn':fildyn,'fildvscf':fildvscf,
                 'outdir':"'%s'"%outdir,'trans':'.True.','ldisp':'.True.','recover':restart,
                 'electron_phonon':"'interpolated'",
                 'nq1':q_mesh_dyn[0],'nq2':q_mesh_dyn[1],'nq3':q_mesh_dyn[2]}
    fs_inputph=make_fstring_obj('inputph',var_inputph,val_inputph,'ph')
    fstream=fstream+fs_inputph
    write_file(fname,fstream)

def make_q2r():
    fname='%s.q2r'%(prefix)
    fstream=''
    fs_input=make_fstring_obj('input',['fildyn','flfrc'],{'fildyn':fildyn,'flfrc':flfrc},'q2r')
    fstream=fstream+fs_input
    write_file(fname,fstream)

def make_matdyn(phband):
    asr='crystal'
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
        fs_klist=k_line_stream(q_mesh_bands,k_list)
        fstream=fstream+fs_klist
    write_file(fname,fstream)

def make_plotband_in(mode,win):
    check_type(mode,bool)
    if mode:
        fname='eband.plotband'
        name=prefix
        fext='band'
        ef=get_ef()
        win_min=win[0]+ef
        win_max=win[1]+ef
    else:
        fname='pband.plotband'
        name='matdyn'
        fext='freq'
        ef=0.
        win_min=0
        win_max=win[1]
    format=(name,fext,win_min,win_max,name,name,ef,nband,ef)
    fstream='%s.%s\n%d %d\n%s.xmgr\n%s.ps\n%f\n%f %f\n'%format
    write_file(fname,fstream)

def main(prefix):
    make_pw_in('scf')
    if sw_bands:
        make_pw_in('bands')
        make_bands_in()
        make_plotband_in(True,eband_win)
        if sw_wan:
            make_pw_in('nscf')
    if sw_ph:
        make_ph_in()
        make_q2r()
        make_matdyn(False)
        make_matdyn(True)
        make_plotband_in(False,pband_win)
if __name__=="__main__":
    main(prefix)
