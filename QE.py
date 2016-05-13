#!/usr/bin/env python
# -*- coding:utf-8
prefix='FeSe'

aa=3.77376; ab=aa; ac=5.52482
zp=0.2652
brav='P'
axis=[aa,ab,ac]
deg=[90,90,90]
atom=['Fe','Se']
atom_position=[[[0.25,0.75,0.],[0.75,0.25,0.]],
               [[0.25,0.25,zp],[0.75,0.75,1-zp]]]

type_xc='pbe'
pot_type=['%s.%s-sp-van','%s.%s-van']
UPF=[pp%(at,type_xc) for pp,at in zip(pot_type,atom)]
psude_dir='/home/Apps/upf_files/'

#===================pw_parameters===============================
k_mesh_scf=10
k_mesh_bands=20
k_mesh_wannier=8
ecut=40.0
conv=1.0d-9
nband=50
deg=0.025
#===================ph_parameters===============================
q_mesh_dyn=[4,4,4]
q_mesh_bands=20
q_mesh_dos=8
#===================Wannier_parameters=============================
nwann=10                     #number of wannier
dis_win=[-3.00,4.00]         #max(min)_window
frz_win=[-0.0,0.0]           #froz_window
projection=[['Fe','d']]      #projections
sw_fs_plot= False            #plot Fermi surface
unk=False                    #Bloch(Wannier)_func

#================physical parameter================================
bohr=0.529177210818; ibohr=1.0/bohr
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
#==================================================================

def write_file(name,stream):
    f=open(name,'r')
    f.write(stream)
    f.close()

def cell_parameter_stream(axis,deg):
    cell_string='cell_parameters\n'
    return cell_string

def k_line_stream(k_points):
    kstring='%d'%len(k_points)
    return k_string

def make_fstring_obj(obj_name,var_list):
    var_name=globals()
    fstring='%s\n'%obj_name
    for vc in var_list:
        fstring=fstring+'%28s = '%vc+str(var_name(vc))+'\n'
    fstring=fstring+'/\n'
    return fstring

def make_pw_in():
    fext='nscf' if calc in ['nscf','bands'] else 'scf'
    fname='%s.%s'%(prefix,fext)
    calculation=calc
    fstream=''
    var_control=['title','calculation','restart_mode','outdir','psude_dir',
                 'prefix','etot_conv_thr','forc_conv_thr','nstep','tstress','tprnfor']
    fs_control=make_fstring_obj('control',var_control)
    fstream=fstream+fs_control
    var_system=[]
    fs_system=make_fstring_obj('system',var_system)
    var_electrons=[]
    fs_electrons=make_fstring_obj('electrons',var_electrons)
    if calc in ['relax','md','vc-relax','vc-md']:
        var_ions=[]
        fs_ions=make_fstring_obj('ions',var_ions)
    if calc in ['vc-relax','vc-md']:
        var_cell=[]
        fs_cell=make_fstring_obj('cell',var_cell)
    fs_atomparam=atomic_parameters_stream()
    if ibrav==0:
        fs_cellparam=cell_parameter_stream()
    write_file(fname,fstream)

def make_bands_in():
    pass

def make_ph_in():
    fname='%s.ph'%(prefix)
    fstream='%s\n'%prefix
    var_inputph=[]
    fs_inputph=make_fstring_obj('inputph',var_sinputph)
    write_file(fname,fstream)

def make_q2r():
    fname='%s.q2r'%(prefix)
    fstream=''
    fildyn='%s.dyn'%prefix
    flfrc='%s.fc'%prefix
    fs_input=make_fstring_obj('input',['fildyn','flfrc'])
    write_file(fname,fstream)

def make_matdyn():
    if phband:
        fext='freq'
        dos=False
        input_list=['flfrc','asr','dos']
    else:
        fext='matdyn'
        dos=True
        nk1=q_mesh_dos
        nk2=q_mesh_dos
        nk3=q_mesh_dos
        input_list=['flfrc','asr','dos','la2F','nk1','nk2','nk3']
    fname='%s.%s'%(prefix,fext)
    fstream=''
    flfrc='%s.fc'%prefix
    asr='crystal'
    fs_input=make_fstring_obj('input',input_list)
    fstream=fstream+fs_input
    if phband:
        fs_klist=k_line_stream()
        fstream=fstream+fs_klist
    write_file(fname,fstream)

def make_plotband_in(mode):
    if mode:
        fname='eband.plotband'
        name=prefix
        fext='band'
        ef=get_ef()
    else:
        fname='pband.plotband'
        name='matdyn'
        fext='freq'
        ef=0.
    format=(name,fext,win_in,win_max,name,name,ef,nband,ef)
    fstream='%s.%s\n%d %d\n%s.xmgr\n%s.ps\n%f\n%f %f\n'%format
    write_file(fname,fstream)
