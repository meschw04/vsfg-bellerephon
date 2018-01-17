import pandas as pd
import numpy as np
import csv
import math
import matplotlib.pyplot as plt

def parse_pdb_file(file_name):
    data = open(file_name,'r').read()
    t = data.split('\n')
    header = []
    body = []
    for i in t:
        if 'ATOM' not in i:
            header.append(i)
        else:
            sub_li = i.split(' ')[0:len(i.split(' '))-3]
            new_li = []
            for j in sub_li:
                if len(j) == 0:
                    pass
                else:
                    try:
                        new_li.append(float(j))
                    except ValueError:
                        new_li.append(j)
            if len(new_li) != 11:
                t = []
                t = t+new_li[:1]
                t.append(new_li[2][:len(new_li[2])-3])
                t.append(new_li[2][len(new_li[2])-3:])
                t = t + new_li[3:]
                new_li = t
            else:
                pass
            body.append(new_li)
    cleaned = []
    for i in body:
        t = []
        for j in i:
            if j != '':
                t.append(j)
            else:
                continue
        if type(t[1]) == float:
            del t[1]
        if t[-1] == 0.0:
            del t[-1]
        if t[-1] == 1.0:
            del t[-1]
        if len(t) == 9:
            t = [t[0],t[1]+t[2],t[3],t[4],t[5],t[6],t[7],t[8]]
        cleaned.append(t)
    body = cleaned
    chain_address = []
    for i in body:
        chain_address.append(i[3])
    chain_address = list(set(chain_address))
    return header, body, chain_address


def position_compare(position_nitrogen,position_carbon):
    '''
    if position_nitrogen == position_carbon:
        C_N_pairs = []
        for i in range(0,len(position_nitrogen)-1):
            coord_m = atom_full_positions[position_nitrogen[i]]
            coord_n = atom_full_positions[position_nitrogen[i+1]]
            new_coord = [coord_m[0]-coord_n[0],\
                         coord_m[1]-coord_n[1],\
                         coord_m[2]-coord_n[2]]
            C_N_pairs.append(new_coord)
        return C_N_pairs
    '''

    C_N_pairs = []
    print(position_nitrogen,position_carbon) #N = oxygen, C = carbon
    if len(position_nitrogen) == len(position_carbon):
        for pos in range(0,len(position_nitrogen)):
            n_atom = atom_full_positions[position_nitrogen[pos]]
            c_atom = atom_full_positions[position_carbon[pos]]
            print(n_atom,position_nitrogen[pos])
            print(c_atom,position_carbon[pos])
            C_N_pairs.append([c_atom[0]-n_atom[0],\
                              c_atom[1]-n_atom[1],\
                              c_atom[2]-n_atom[2]])
    else:
        for i in range(0,len(max(position_nitrogen,position_carbon))):
            if abs(position_nitrogen[i]-position_carbon[i]) < 3:
                pass
            else:
                if len(position_carbon) < len(position_nitrogen):
                    del position_nitrogen[i]
                else:
                    del position_carbon[i]
        if len(position_nitrogen) != len(position_carbon):
            if len(position_nitrogen) < len(position_carbon):
                del position_carbon[-1]
            else:
                del position_nitrogen[-1]
        for pos in range(0,len(position_nitrogen)):
            n_atom = atom_full_positions[position_nitrogen[pos]]
            c_atom = atom_full_positions[position_carbon[pos]]
            C_N_pairs.append([c_atom[0]-n_atom[0],\
                              c_atom[1]-n_atom[1],\
                              c_atom[2]-n_atom[2]])
    return C_N_pairs


#header, body, chain_address = parse_pdb_file('/Users/mschwart/vsfg-bellerephon/sfg2/output_CHAIN_A.pdb') #PRepped PDB

header, body, chain_address = parse_pdb_file('/Users/mschwart/vsfg-bellerephon/test-16jan/output.pdb') #PRepped PDB
#print(body)
#header, body, chain_address = parse_pdb_file('C:\Users\Marcus\Documents\Fall 2017\sfg2\output.pdb')

#header, body, chain_address = parse_pdb_file('C:\Users\Marcus\Documents\Fall 2017\sfg2\output_TZ.pdb')


#This is the cR from the Mathematica code
trans_dip_mom = pd.read_csv('/Users/mschwart/vsfg-bellerephon/test-16jan/output_prepped_transdipmom_Hamm.txt',header=None) #Trans Dipole Moment
trans = []
for i in trans_dip_mom[0].tolist():
    t = i.split(' ')
    l = []
    for k in t:
        if len(k) == 0:
            pass
        else:
            l.append(-float(k))
    trans.append(l)


#EVAL is for the eigenvalues
eigenvalues = pd.read_csv('/Users/mschwart/vsfg-bellerephon/test-16jan/output_prepped_eval.txt',header=None) #Eigenvalues
eig_val = []
for i in eigenvalues[0].tolist():
    eig_val.append(float(i)+1660)
eigenvalues = eig_val

#EVEC is for the eigenvectors
eigenvectors = pd.read_csv('/Users/mschwart/vsfg-bellerephon/test-16jan/output_prepped_evec.txt',header=None) #Eigenvectors

full_eig = []
for i in eigenvectors[0].tolist():
    t=[]
    for num in i.split(' '):
        if len(num) == 0:
            pass
        else:
            t.append(float(num))
    full_eig.append(t)
eigenvectors = full_eig #This is a 288 * 288 matrix.

atom_full_positions = []
for i in body:
    atom_full_positions.append([float(i[5]),float(i[6]),float(i[7])])
#    atom_full_positions.append([float(i[5]),float(i[6]),float(i[7])])



def atom_analysis(header,body):
    atom_chars = [i[1] for i in body]
    position_nitrogen = []
    position_carbon = []
    position_oxygen = []
    position_hydrogen = []
    position_carbon_alpha = []
    
    for atom in range(0,len(atom_chars)):
        if atom_chars[atom] == 'H':
            position_hydrogen.append(atom)
        elif atom_chars[atom] == 'O':
            position_oxygen.append(atom)
    
    for atom in range(0,len(atom_chars)):
        if len(position_hydrogen) != len(position_nitrogen):
            if atom_chars[atom] == 'N' and \
            abs(position_hydrogen[len(position_nitrogen)]-atom)<3:
                position_nitrogen.append(atom)
        if len(position_oxygen) != len(position_carbon):
            if atom_chars[atom] == 'C' and \
            abs(position_oxygen[len(position_carbon)]-atom)<3:
                position_carbon.append(atom)
        if len(position_hydrogen) != len(position_carbon_alpha):
            if atom_chars[atom] == 'CA' and \
            abs(position_hydrogen[len(position_carbon_alpha)]-atom)<3:
                position_carbon_alpha.append(atom)
    
        else:
            pass
    print(position_oxygen,'\n',position_carbon)
    C_N_pairs = position_compare(position_nitrogen,position_carbon)
    C_O_pairs = position_compare(position_oxygen,position_carbon)
    print(C_O_pairs)
    c = []
    for i in C_O_pairs:
        norm_int = (i[0]**2 + i[1]**2 + i[2]**2)**0.5
        new_norm = []
        for k in i:
            try:
                new_norm.append(k/norm_int)
            except ZeroDivisionError:
                new_norm.append(0.0)
        c.append(new_norm)
    print(c)
    a = []
    for i in range(0,len(c)):
        numer = list(np.cross(c[i],C_N_pairs[i]))
        denom = (numer[0]**2+numer[1]**2+numer[2]**2)**0.5
        inter = []
        for j in numer:
            try:
                inter.append(j/denom)
            except ZeroDivisionError:
                inter.append(0.0)
        a.append(inter)
    
    
    b = []
    for i in range(0,len(c)):
        crossed = list(np.cross(c[i],a[i]))
        b.append(crossed)
    print(b)
    Raman_tensor_1_mode = [[0.05,0.0,0.0],[0.0,0.2,0.0],[0.0,0.0,1.0]]
    
    alpha_Z = []
    for i in range(0,len(b)):
        trig_cos = math.cos(146*math.pi/180)
        trig_sin = math.sin(146*math.pi/180)
        man_app = []
        for j in range(0,len(c[i])):
            man_app.append((trig_cos*c[i][j])-(trig_sin*b[i][j]))
        alpha_Z.append(man_app)
    
    alpha_X = []
    for i in range(0,len(alpha_Z)):
        numer = list(np.cross(alpha_Z[i],C_N_pairs[i]))
        denom = (numer[0]**2+numer[1]**2+numer[2]**2)**0.5
        inter = []
        for j in numer:
            inter.append(j/denom)
        alpha_X.append(inter)
    alpha_Y = []
    for i in range(0,len(alpha_X)):
        crossed = list(np.cross(alpha_Z[i],alpha_X[i]))
        alpha_Y.append(crossed)
    
    norm_alpha_X = []
    for i in alpha_X:
        norm_alpha_X.append((i[0]**2+i[1]**2+i[2]**2)**0.5)
    norm_alpha_Y = []
    for i in alpha_Y:
        norm_alpha_Y.append((i[0]**2+i[1]**2+i[2]**2)**0.5)
    norm_alpha_Z = []
    for i in alpha_Z:
        norm_alpha_Z.append((i[0]**2+i[1]**2+i[2]**2)**0.5)
    
    Raman_Tensor_Table_AA = []
    for i in range(0,len(alpha_X)):
        t = (((alpha_X[i][0]/norm_alpha_X[i])**2)*Raman_tensor_1_mode[0][0])+\
            (((alpha_Y[i][0]/norm_alpha_Y[i])**2)*Raman_tensor_1_mode[1][1])+\
            (((alpha_Z[i][0]/norm_alpha_Z[i])**2)*Raman_tensor_1_mode[2][2])
        Raman_Tensor_Table_AA.append(t)
    
    Raman_Tensor_Table_BB = []
    for i in range(0,len(alpha_X)):
        t = (((alpha_X[i][1]/norm_alpha_X[i])**2)*Raman_tensor_1_mode[0][0])+\
            (((alpha_Y[i][1]/norm_alpha_Y[i])**2)*Raman_tensor_1_mode[1][1])+\
            (((alpha_Z[i][1]/norm_alpha_Z[i])**2)*Raman_tensor_1_mode[2][2])
        Raman_Tensor_Table_BB.append(t)
    
    Raman_Tensor_Table_CC = []
    for i in range(0,len(alpha_X)):
        t = (((alpha_X[i][2]/norm_alpha_X[i])**2)*Raman_tensor_1_mode[0][0])+\
            (((alpha_Y[i][2]/norm_alpha_Y[i])**2)*Raman_tensor_1_mode[1][1])+\
            (((alpha_Z[i][2]/norm_alpha_Z[i])**2)*Raman_tensor_1_mode[2][2])
        Raman_Tensor_Table_CC.append(t)
    
    Raman_Tensor_Table_BC = []
    for i in range(0,len(alpha_X)):
        t = (((alpha_X[i][1]/norm_alpha_X[i])*(alpha_X[i][2]/norm_alpha_X[i]))*Raman_tensor_1_mode[0][0])+\
            (((alpha_Y[i][1]/norm_alpha_Y[i])*(alpha_Y[i][2]/norm_alpha_Y[i]))*Raman_tensor_1_mode[1][1])+\
            (((alpha_Z[i][1]/norm_alpha_Z[i])*(alpha_Z[i][2]/norm_alpha_Z[i]))*Raman_tensor_1_mode[2][2])
        Raman_Tensor_Table_BC.append(t)
    
    Raman_Tensor_Table_AB = []
    for i in range(0,len(alpha_X)):
        t = (((alpha_X[i][0]/norm_alpha_X[i])*(alpha_X[i][1]/norm_alpha_X[i]))*Raman_tensor_1_mode[0][0])+\
            (((alpha_Y[i][0]/norm_alpha_Y[i])*(alpha_Y[i][1]/norm_alpha_Y[i]))*Raman_tensor_1_mode[1][1])+\
            (((alpha_Z[i][0]/norm_alpha_Z[i])*(alpha_Z[i][1]/norm_alpha_Z[i]))*Raman_tensor_1_mode[2][2])
        Raman_Tensor_Table_AB.append(t)
    
    Raman_Tensor_Table_AC = []
    for i in range(0,len(alpha_X)):
        t = (((alpha_X[i][0]/norm_alpha_X[i])*(alpha_X[i][2]/norm_alpha_X[i]))*Raman_tensor_1_mode[0][0])+\
            (((alpha_Y[i][0]/norm_alpha_Y[i])*(alpha_Y[i][2]/norm_alpha_Y[i]))*Raman_tensor_1_mode[1][1])+\
            (((alpha_Z[i][0]/norm_alpha_Z[i])*(alpha_Z[i][2]/norm_alpha_Z[i]))*Raman_tensor_1_mode[2][2])
        Raman_Tensor_Table_AC.append(t)
    
    Raman_Tensor_Table = []
    for i in range(0,len(Raman_Tensor_Table_AA)):
        Raman_Tensor_Table.append([[Raman_Tensor_Table_AA[i],Raman_Tensor_Table_AB[i],Raman_Tensor_Table_AC[i]],\
                                   [Raman_Tensor_Table_AB[i],Raman_Tensor_Table_BB[i],Raman_Tensor_Table_BC[i]],\
                                   [Raman_Tensor_Table_AC[i],Raman_Tensor_Table_BC[i],Raman_Tensor_Table_CC[i]]])
    
    Iir1mode = []
    for k in range(0,len(eigenvectors)):
        t1 = 0
        t2 = 0
        t3 = 0
        for i in range(0,len(eigenvectors[k])):
            q = eigenvectors[k][i]
            li = trans[i]
            t1 = t1+q*li[0]
            t2 = t2+q*li[1]
            t3 = t3+q*li[2]
        Iir1mode.append([t1,t2,t3])
    
    
    Iir1modeAbs = []
    for k in range(0,len(eigenvectors)): ###i!!!
        fff = []
        for i in range(0,len(eigenvectors[k])): ###j!!!
            li = trans[i]
            q = eigenvectors[k][i]
            fff.append([li[0]*q,li[1]*q,li[2]*q])
        fff = [sum(i) for i in zip(*fff)]
        for i in fff:
            Iir1modeAbs.append(fff[0]**2+fff[1]**2+fff[2]**2)
    t = []
    for i in range(0,len(Iir1modeAbs)):
        if i%3 == 0:
            t.append(Iir1modeAbs[i])
        else:
            pass
    Iir1modeAbs=t
    return Iir1mode, Iir1modeAbs, Raman_Tensor_Table
    


Iir1mode, Iir1modeAbs, Raman_Tensor_Table = atom_analysis(header,body)


##CONSTANTS##
Gamma_IR_Sticks=0.01
Gamma_IR=10.0
Omega_offset_IR=0.0



def ir_stick_spectrum_compute(omega,eigenvalues,Iir1modeAbs):
    q = []
    for i in range(0,len(Iir1modeAbs)):
        q.append((Iir1modeAbs[i]*Gamma_IR_Sticks/np.pi())/((omega - eigenvalues[i]-Omega_offset_IR)**2+Gamma_IR_Sticks**2))
    IR_Stick_Spectrum_Omega = sum(q)
    return IR_Stick_Spectrum_Omega

def ir_spectrum_compute(omega,eigenvalues,Iir1modeAbs):
    q = []
    for i in range(0,len(eigenvalues)):

        q.append((Iir1modeAbs[i]*Gamma_IR/float(math.pi))/((omega - eigenvalues[i]-Omega_offset_IR)**2+Gamma_IR**2))
    IR_Spectrum_Omega = sum(q)
    return IR_Spectrum_Omega

IRspec = []
for i in range(3100,3501):
    IRspec.append(ir_spectrum_compute(0.5*float(i),eigenvalues,Iir1modeAbs))
#print(IRspec)

OmegaVal=[x * 0.5 for x in range(3100, 3501)]

plt.figure(figsize=(8,8))
plt.title('IR Curve')
plt.plot(OmegaVal,IRspec)
plt.xlabel('Wavelength (nm)')
plt.ylabel('Intensity (normalized)')
plt.show()

Iraman1mode = []
Iraman1modeIsotropic= []
for k in range(0,len(eigenvectors)): #row of eigenvector matrix
    a=0.0
    b=0.0
    c=0.0
    d=0.0
    e=0.0
    f=0.0
    g=0.0
    h=0.0
    i=0.0
    for v in range(0,len(eigenvectors[k])): #value within eigenvector row
        i_j = eigenvectors[k][v]
        Raman_matrix = Raman_Tensor_Table[v] #3x3 matrix
        a = a + float(Raman_matrix[0][0]*i_j)
        b = b + float(Raman_matrix[0][1]*i_j)
        c = b + float(Raman_matrix[0][2]*i_j)
        d = b + float(Raman_matrix[1][0]*i_j)
        e = b + float(Raman_matrix[1][1]*i_j)
        f = b + float(Raman_matrix[1][2]*i_j)
        g = b + float(Raman_matrix[2][0]*i_j)
        h = b + float(Raman_matrix[2][1]*i_j)
        i = b + float(Raman_matrix[2][2]*i_j)
    Iraman1mode.append([[a,b,c],[d,e,f],[g,h,i]])
    Iraman1modeIsotropic.append(((a*e*i)**2)/9.0) #a*e*i = Trace of the matrix

Iraman1modeAnisotropic = [] #Maybe want to incorporate this into the above code block later...
for i in Iraman1mode:
    tmp1 = .5*(i[0][0]-i[1][1])**2
    tmp2 = .5*(i[1][1]-i[2][2])**2
    tmp3 = .5*(i[2][2]-i[0][0])**2
    tmp3 = .5*(i[2][2]-i[0][0])**2
    tmp4 = 0.75*(i[0][1]-i[1][0])**2
    tmp5 = 0.75*(i[1][2]-i[2][1])**2
    tmp6 = 0.75*(i[0][2]-i[2][0])**2
    Iraman1modeAnisotropic.append(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6)

Omega_offset_IR = 0.0


def Raman_Spectrum_Isotropic(omega,Gamma_Raman):
    term = 0.0
    for i in range(0,len(Iraman1modeIsotropic)):
        term = term+(Iraman1modeIsotropic[i]*Gamma_Raman/math.pi)/(omega-eigenvalues[i]-Omega_offset_IR)**2
    return term

def Raman_Spectrum_Anisotropic(omega,Gamma_Raman):
    term = 0.0
    for i in range(0,len(Iraman1modeIsotropic)):
        term = term+(Iraman1modeAnisotropic[i]*Gamma_Raman/math.pi)/(omega-eigenvalues[i]-Omega_offset_IR)**2
    return term
    


#CONSTANT!!!#
Capital_Gamma_Raman = 9.0



linspace = range(1575,1751)
t = []
s = []
for i in linspace:
    t.append(Raman_Spectrum_Isotropic(float(i),Capital_Gamma_Raman))
    s.append(Raman_Spectrum_Anisotropic(float(i),Capital_Gamma_Raman))
plt.figure(figsize=(8,8))
plt.title('Raman Spectra')
plt.plot(linspace,t)
plt.xlabel('Wave Number (cm-1)')
plt.ylabel('Intensity (normalized)')
plt.show()

plt.figure(figsize=(8,8))
plt.title('Raman Spectra (Anisotropic)')
plt.plot(linspace,s)
plt.xlabel('Wave Number (cm-1)')
plt.ylabel('Intensity (normalized)')
plt.show()

Isfg1mode = []
for i in range(0,len(Iraman1mode)):
    t = Iraman1mode[i]
    s = Iir1mode[i]
    build = [[[t[0][0]*s[0],t[0][0]*s[1],t[0][0]*s[2]],
               [t[0][1]*s[0],t[0][1]*s[1],t[0][1]*s[2]],
               [t[0][0]*s[0],t[0][0]*s[1],t[0][0]*s[2]]],
               [[t[0][0]*s[0],t[0][0]*s[1],t[0][0]*s[2]],
               [t[0][0]*s[0],t[0][0]*s[1],t[0][0]*s[2]],
               [t[0][0]*s[0],t[0][0]*s[1],t[0][0]*s[2]]],
               [[t[0][0]*s[0],t[0][0]*s[1],t[0][0]*s[2]],
               [t[0][0]*s[0],t[0][0]*s[1],t[0][0]*s[2]],
               [t[0][0]*s[0],t[0][0]*s[1],t[0][0]*s[2]]]]
    Isfg1mode.append(build)

Chi2XXZ1Molecule = [x[0][0][2] for x in Isfg1mode]

Chi2YYZ1Molecule = [x[1][1][2] for x in Isfg1mode]

Chi2ZZZ1Molecule = [x[2][2][2] for x in Isfg1mode]

Chi2XZX1Molecule = [x[0][2][0] for x in Isfg1mode]

Chi2YZY1Molecule = [x[1][2][1] for x in Isfg1mode]

Chi2ZZX1Molecule = [x[2][2][0] for x in Isfg1mode]

Chi2ZZY1Molecule = [x[2][2][1] for x in Isfg1mode]

xxz = [x[0][0][2] for x in Isfg1mode]

yyz = [x[1][1][2] for x in Isfg1mode]

zzz = [x[2][2][2] for x in Isfg1mode]

xzx = [x[0][2][0] for x in Isfg1mode]

yzy = [x[1][2][1] for x in Isfg1mode]

zzx = [x[2][2][0] for x in Isfg1mode]

zzy = [x[2][2][1] for x in Isfg1mode]

yzz = [x[1][2][2] for x in Isfg1mode]

zyz = [x[2][1][2] for x in Isfg1mode]

xzz = [x[0][2][2] for x in Isfg1mode]

zxz = [x[2][0][2] for x in Isfg1mode]

zxx = [x[2][0][0] for x in Isfg1mode]

zyy = [x[2][1][1] for x in Isfg1mode]

xyz = [x[0][1][2] for x in Isfg1mode]

xzy = [x[0][2][1] for x in Isfg1mode]

yxz = [x[1][0][2] for x in Isfg1mode]

yzx = [x[1][2][0] for x in Isfg1mode]

zxy = [x[2][0][1] for x in Isfg1mode]

zyx = [x[2][1][0] for x in Isfg1mode]

xxy = [x[0][0][1] for x in Isfg1mode]

xyx = [x[0][1][0] for x in Isfg1mode]

yxx = [x[1][0][0] for x in Isfg1mode]

yzz = [x[1][2][2] for x in Isfg1mode]

zyz = [x[2][1][2] for x in Isfg1mode]

xyy = [x[0][1][1] for x in Isfg1mode]

yxy = [x[1][0][1] for x in Isfg1mode]

yyx = [x[1][1][0] for x in Isfg1mode]

yyy = [x[1][1][1] for x in Isfg1mode]

xxx = [x[0][0][0] for x in Isfg1mode]

xyy = [x[0][1][1] for x in Isfg1mode]

yxy = [x[1][0][1] for x in Isfg1mode]

yyx = [x[1][1][0] for x in Isfg1mode]


Ns = 1.0

def Chi2DeltaDistZZZEnsembleDeltaDist(theta_center): #theta_center must be in degrees, not radians!
    rad_t = theta_center*math.pi/180
    returned_val = list(np.zeros(len(zzz)))
    mark1 = math.cos(rad_t)**3
    returned_val = [x + mark1*y for x, y in zip(returned_val, zzz)]
    mark2 = 0.5*(math.sin(rad_t)**2)*math.cos(rad_t)
    returned_val = [x + mark2*(y+z+w) for x, y, z, w in zip(returned_val, yyz, yzy, zyy)]
    mark3 = 0.5*(math.sin(rad_t)**2)*math.cos(rad_t)
    returned_val = [x + mark3*(y+z+w) for x, y, z, w in zip(returned_val, xxz, xzx, zxx)]
    returned_val = [Ns*i for i in returned_val]
    return returned_val

def Chi2DeltaDistXZXEnsembleDeltaDist(theta_center):
    rad_t = theta_center*math.pi/180
    returned_val = list(np.zeros(len(zzz)))
    mark1 = (math.sin(rad_t)**2)*math.cos(rad_t)
    returned_val = [x + mark1*y for x, y in zip(returned_val, zzz)]
    mark2 = math.cos(rad_t)
    returned_val = [x + mark2*(y+z) for x, y, z in zip(returned_val, xzx, yzy)]
    mark3 = 0.5*math.cos(rad_t)*math.sin(rad_t)**2
    returned_val = [x + mark3*(y+z+w) for x, y, z, w in zip(returned_val, yyz, yzy, zyy)]
    mark4 = -0.5*math.cos(rad_t)*math.sin(rad_t)**2
    returned_val = [x + mark4*(y+z+w) for x, y, z, w in zip(returned_val, xxz, xzx, zxx)]
    returned_val = [0.5*Ns*i for i in returned_val]
    return returned_val

def Chi2DeltaDistXXZEnsembleDeltaDist(theta_center):
    rad_t = theta_center*math.pi/180
    returned_val = list(np.zeros(len(zzz)))
    mark1 = (math.sin(rad_t)**2)*math.cos(rad_t)
    returned_val = [x + mark1*y for x, y in zip(returned_val, zzz)]
    mark2 = math.cos(rad_t)
    returned_val = [x + mark2*(y+z) for x, y, z in zip(returned_val, xxz, yyz)]
    mark3 = 0.5*(math.sin(rad_t)**2)*math.cos(rad_t)
    returned_val = [x + mark3*(y+z+w) for x,y,z,w in zip(returned_val, yyz, yzy, zyy)]
    mark4 = -0.5*(math.sin(rad_t)**2)*math.cos(rad_t)
    returned_val = [x + mark4*(y+z+w) for x,y,z,w in zip(returned_val, xxz, xzx, zxx)]
    returned_val = [0.5*Ns*i for i in returned_val]
    return returned_val

def Chi2DeltaDistZXXEnsembleDeltaDist(theta_center):
    rad_t = theta_center*math.pi/180
    returned_val = list(np.zeros(len(zzz)))
    mark1 = (math.sin(rad_t)**2)*math.cos(rad_t)
    returned_val = [x + mark1*y for x, y in zip(returned_val, zzz)]
    mark2 = math.cos(rad_t)
    returned_val = [x + mark2*(y+z) for x, y, z in zip(returned_val, zxx, zyy)]
    mark3 = 0.5*(math.sin(rad_t)**2)*math.cos(rad_t)
    returned_val = [x + mark3*(y+z+w) for x,y,z,w in zip(returned_val, yyz, yzy, zyy)]
    mark4 = -0.5*(math.sin(rad_t)**2)*math.cos(rad_t)
    returned_val = [x + mark4*(y+z+w) for x,y,z,w in zip(returned_val, xxz, xzx, zxx)]
    returned_val = [0.5*Ns*i for i in returned_val]
    return returned_val



Gamma_Exc = 7.5




#ISSUE HERE!
def Chi2DeltaDistXXZ(omega, theta_center, omega_offset, capital_gamma):
    #omega is float, theta_center is list of same length as eigenvalues, 
    #omega_offset is float, capital_gamma is float.
    re_val = 0
    im_val = 0
    for j in range(0,len(eigenvalues)):
        denom = ((eigenvalues[j]+omega_offset-omega)**2)+(Gamma_Exc+capital_gamma)**2
        re_num = (eigenvalues[j]+omega_offset-omega)*((Gamma_Exc+capital_gamma)**0.5)*\
                 Chi2DeltaDistXXZEnsembleDeltaDist(theta_center)[j]
        re_val = re_val+(re_num/denom)
        im_num = ((Gamma_Exc+capital_gamma)**1.5)*\
                 Chi2DeltaDistXXZEnsembleDeltaDist(theta_center)[j]
        im_val = im_val+(im_num/denom)
    return (re_val,im_val)
    
def Chi2DeltaDistZZZ(omega, theta_center, omega_offset, capital_gamma):
    #omega is float, theta_center is list of same length as eigenvalues, 
    #omega_offset is float, capital_gamma is float.
    re_val = 0
    im_val = 0
    for j in range(0,len(eigenvalues)):
        denom = ((eigenvalues[j]+omega_offset-omega)**2)+(Gamma_Exc+capital_gamma)**2
        re_num = (eigenvalues[j]+omega_offset-omega)*((Gamma_Exc+capital_gamma)**0.5)*\
                 Chi2DeltaDistZZZEnsembleDeltaDist(theta_center)[j]
        re_val = re_val+(re_num/denom)
        im_num = ((Gamma_Exc+capital_gamma)**1.5)*\
                 Chi2DeltaDistZZZEnsembleDeltaDist(theta_center)[j]
        im_val = im_val+(im_num/denom)
    return (re_val,im_val)

def Chi2DeltaDistXZX(omega, theta_center, omega_offset, capital_gamma):
    #omega is float, theta_center is list of same length as eigenvalues, 
    #omega_offset is float, capital_gamma is float.
    re_val = 0
    im_val = 0
    for j in range(0,len(eigenvalues)):
        denom = ((eigenvalues[j]+omega_offset-omega)**2)+(Gamma_Exc+capital_gamma)**2
        re_num = (eigenvalues[j]+omega_offset-omega)*((Gamma_Exc+capital_gamma)**0.5)*\
                 Chi2DeltaDistXZXEnsembleDeltaDist(theta_center)[j]
        re_val = re_val+(re_num/denom)
        im_num = ((Gamma_Exc+capital_gamma)**1.5)*\
                 Chi2DeltaDistXZXEnsembleDeltaDist(theta_center)[j]
        im_val = im_val+(im_num/denom)
    return (re_val,im_val)
    
def Chi2DeltaDistZXX(omega, theta_center, omega_offset, capital_gamma):
    #omega is float, theta_center is list of same length as eigenvalues, 
    #omega_offset is float, capital_gamma is float.
    re_val = 0
    im_val = 0
    for j in range(0,len(eigenvalues)):
        denom = ((eigenvalues[j]+omega_offset-omega)**2)+(Gamma_Exc+capital_gamma)**2
        re_num = (eigenvalues[j]+omega_offset-omega)*((Gamma_Exc+capital_gamma)**0.5)*\
                 Chi2DeltaDistZXXEnsembleDeltaDist(theta_center)[j]
        re_val = re_val+(re_num/denom)
        im_num = ((Gamma_Exc+capital_gamma)**1.5)*\
                 Chi2DeltaDistZXXEnsembleDeltaDist(theta_center)[j]
        im_val = im_val+(im_num/denom)
    return (re_val,im_val)
'''
def Chi2DeltaDistZXZ(omega, theta_center, omega_offset, capital_gamma):
    #omega is float, theta_center is list of same length as eigenvalues, 
    #omega_offset is float, capital_gamma is float.
    re_val = 0
    im_val = 0
    for j in range(0,len(eigenvalues)):
        denom = ((eigenvalues[j]+omega_offset-omega)**2)+(Gamma_Exc+capital_gamma)**2
        re_num = (eigenvalues[j]+omega_offset-omega)*((Gamma_Exc+capital_gamma)**0.5)*\
                 Chi2DeltaDistZXZEnsembleDeltaDist(theta_center)[j]
        re_val = re_val+(re_num/denom)
        im_num = ((Gamma_Exc+capital_gamma)**1.5)*\
                 Chi2DeltaDistZXZEnsembleDeltaDist(theta_center)[j]
        im_val = im_val+(im_num/denom)
    return (re_val,im_val)
'''

incIR=40*math.pi/180.0
incVIS=36.*math.pi/180.0
k1 = 2.*math.pi/800.0
k2 = 2.*math.pi/6000.0
n1SF=1.0
n2SF=1.331
n1VIS=1.0
n2VIS=1.329
n1IR=1.0
n2IR=1.315
Theta_Center=0.1
#Capital_Gamma=0.1
Capital_Gamma = 10.0
Omega_offset=0.0
Omega_CO=1728.0
Capital_Gamma_CO=0.0
fCO = 0.0
ScaleFactorXXZ=0.043578
Anr=0
Phi_nr=0.0
npSF = 1.18
npSFCO = 1.331
npIR = 1.18
npIRCO = 1.315


def Chi2DeltaDistXXZ_NR(omega, theta_center, omega_offset, capital_gamma,A_nr,Phi_nr):
    (re,im) = Chi2DeltaDistXXZ(omega,theta_center,omega_offset,capital_gamma)
    re = re-A_nr*math.exp(Phi_nr/math.pi)
    return (re,im)
    
def Chi2DeltaDistZZZ_NR(omega, theta_center, omega_offset, capital_gamma,A_nr,Phi_nr):
    (re,im) = Chi2DeltaDistZZZ(omega,theta_center,omega_offset,capital_gamma)
    re = re-A_nr*math.exp(Phi_nr/math.pi)
    return (re,im)

def Chi2DeltaDistXZX_NR(omega, theta_center, omega_offset, capital_gamma,A_nr,Phi_nr):
    (re,im) = Chi2DeltaDistXZX(omega,theta_center,omega_offset,capital_gamma)
    re = re-A_nr*math.exp(Phi_nr/math.pi)
    return (re,im)

def Chi2DeltaDistZXX_NR(omega, theta_center, omega_offset, capital_gamma,A_nr,Phi_nr):
    (re,im) = Chi2DeltaDistZXX(omega,theta_center,omega_offset,capital_gamma)
    re = re-A_nr*math.exp(Phi_nr/math.pi)
    return (re,im)

def Chi2DeltaDistZXZ_NR(omega, theta_center, omega_offset, capital_gamma,A_nr,Phi_nr):
    (re,im) = Chi2DeltaDistZXZ(omega,theta_center,omega_offset,capital_gamma)
    re = re-A_nr*math.exp(Phi_nr/math.pi)
    return (re,im)





def Chi2DeltaDistXXZ_NR_LabFrameTotal(omega, theta_center, omega_offset, capital_gamma,A_nr,\
                                      Phi_nr, scale_factor_XXZ):
    (re,im) = Chi2DeltaDistXXZ(omega,theta_center,omega_offset,capital_gamma)
    re = (scale_factor_XXZ*re)-A_nr*math.exp(Phi_nr/math.pi)
    im = (scale_factor_XXZ*im)
    return (re,im)

def Chi2DeltaDistZZZ_NR_LabFrameTotal(omega, theta_center, omega_offset, capital_gamma,A_nr,\
                                      Phi_nr, scale_factor_ZZZ):
    (re,im) = Chi2DeltaDistZZZ(omega,theta_center,omega_offset,capital_gamma)
    re = (scale_factor_ZZZ*re)-A_nr*math.exp(Phi_nr/math.pi)
    im = (scale_factor_ZZZ*im)
    return (re,im)

def Chi2DeltaDistXZX_NR_LabFrameTotal(omega, theta_center, omega_offset, capital_gamma,A_nr,\
                                      Phi_nr, scale_factor_XZX):
    (re,im) = Chi2DeltaDistXZX(omega,theta_center,omega_offset,capital_gamma)
    re = (scale_factor_XZX*re)-A_nr*math.exp(Phi_nr/math.pi)
    im = (scale_factor_XZX*im)
    return (re,im)

def Chi2DeltaDistZXX_NR_LabFrameTotal(omega, theta_center, omega_offset, capital_gamma,A_nr,\
                                      Phi_nr, scale_factor_ZXX):
    (re,im) = Chi2DeltaDistZXX(omega,theta_center,omega_offset,capital_gamma)
    re = (scale_factor_ZXX*re)-A_nr*math.exp(Phi_nr/math.pi)
    im = (scale_factor_ZXX*im)
    return (re,im)

def Chi2DeltaDistZXZ_NR_LabFrameTotal(omega, theta_center, omega_offset, capital_gamma,A_nr,\
                                      Phi_nr, scale_factor_ZXZ):
    (re,im) = Chi2DeltaDistZXZ(omega,theta_center,omega_offset,capital_gamma)
    re = (scale_factor_ZXZ*re)-A_nr*math.exp(Phi_nr/math.pi)
    im = (scale_factor_ZXZ*im)
    return (re,im)

def ImChi2XXZ_NR_LabFrameTotal(omega,theta_center,omega_offset,capital_gamma, A_nr,\
                               Phi_nr, scale_factor_XXZ):
    (re,im) = Chi2DeltaDistXXZ_NR_LabFrameTotal(omega, theta_center, omega_offset,\
                                                capital_gamma,A_nr,\
                                                Phi_nr, scale_factor_XXZ)
    return im

def ReChi2XXZ_NR_LabFrameTotal(omega,theta_center,omega_offset,capital_gamma, A_nr,\
                               Phi_nr, scale_factor_XXZ):
    (re,im) = Chi2DeltaDistXXZ_NR_LabFrameTotal(omega, theta_center, omega_offset,\
                                                capital_gamma,A_nr,\
                                                Phi_nr, scale_factor_XXZ)
    return re


incSF = np.arcsin((k1*math.sin(incVIS)+k2*math.sin(incIR))/(k1+k2))

def Lxx(n1,n2,inc,ref):
    return (2*n1*math.cos(ref))/(n1*math.cos(ref) + n2*math.cos(inc))
def Lyy(n1,n2,inc,ref):
    return (2*n1*math.cos(inc))/(n1*math.cos(inc) + n2*math.cos(ref))

def LzzSF(n1,n2,npSF, inc, ref):
    return (2*n2*math.cos(inc))/(n1*math.cos(ref) + n2*math.cos(inc)) * (n1/npSF)**2
def LzzVIS(n1,n2,npVIS, inc, ref):
    return (2*n2*math.cos(inc))/(n1*math.cos(ref) + n2*math.cos(inc)) * (n1/npVIS)**2
def LzzIR(n1,n2,npIR, inc, ref):
    return (2*n2*math.cos(inc))/(n1*math.cos(inc) + n2*math.cos(inc)) * (n1/npIR)**2

def Lx(n1,n2,inc):
    return Lxx(n1,n2,inc,np.arcsin(n1/n2*math.sin(inc)))
def Ly(n1,n2,inc):
    return Lyy(n1,n2,inc,np.arcsin(n1/n2*math.sin(inc)))

def LzSF(n1,n2,npSF,inc):
    return LzzSF(n1,n2,npSF, inc,np.arcsin(n1/n2*math.sin(inc)))
def LzVIS(n1,n2,npVIS,inc):
    return LzzVIS(n1,n2,npVIS, inc,np.arcsin(n1/n2*math.sin(inc)))
def LzIR(n1,n2,npIR,inc):
    return LzzIR(n1,n2,npIR, inc,np.arcsin(n1/n2*math.sin(inc)))

def Lpppp(npSF, npVIS, npIR):
    return math.sin(incIR)*math.sin(incVIS)*math.sin(incSF)*\
            LzSF(n1SF,n2SF,npSF,incSF)*LzVIS(n1VIS,n2VIS,npVIS,incVIS)*\
                LzIR(n1IR,n2IR,npIR,incIR)

def Lpssp(npIR):
    return math.sin(incIR)*math.cos(incVIS)*math.cos(incSF)*\
           Lx(n1SF,n2SF,incSF)*Lx(n1VIS,n2VIS,incVIS)*LzIR(n1IR,n2IR,npIR,incIR)

def Lpsps(npVIS):
    return -1*math.cos(incIR)*math.sin(incVIS)*math.cos(incSF)*\
            Lx(n1SF,n2SF,incSF)*LzVIS(n1VIS,n2VIS,npVIS,incVIS)\
              *Lx(n1IR,n2IR,incIR)

def Lppss(npSF):
    return math.cos(incIR)*math.cos(incVIS)*math.sin(incSF)*\
            LzSF(n1SF,n2SF,npSF,incSF)*Lx(n1VIS,n2VIS,incVIS)*\
                Lx(n1IR,n2IR,incIR)

def Lssp(npIR):
    return math.sin(incIR)*Ly(n1SF,n2SF, incSF)*\
        Ly(n1VIS,n2VIS,incVIS)*LzIR(n1IR,n2IR,npIR,incIR)

def Lsps(npVIS):
    return math.sin(incVIS)*Ly(n1SF,n2SF, incSF)*\
        LzVIS(n1VIS,n2VIS,npVIS,incVIS)*Ly(n1IR,n2IR,incIR)

def Lpss(npSF):
    return math.sin(incSF)*LzSF(n1SF,n2SF,npSF,incSF)*Ly(n1VIS,n2VIS,incVIS)*Ly(n1IR,n2IR,incIR)

def Chi2DeltaDist_ppp(omega,theta_center,omega_offset,capital_gamma,npSF,npVIS,npIR,Anr,Phi_nr):
    #TAKING ONLY THE REAL PART FOR NOW!
    return Chi2DeltaDistXXZ_NR(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr)[0]*Lpssp(npIR) + \
        Chi2DeltaDistXZX_NR(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr)[0]*Lpsps(npVIS) + \
        Chi2DeltaDistZXX_NR(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr)[0]*Lppss(npSF) + \
        Chi2DeltaDistZZZ_NR(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr)[0]*Lpppp(npSF,npVIS,npIR)
                           
def Chi2DeltaDist_ssp(omega,theta_center,omega_offset,capital_gamma,npIR,Anr,Phi_nr):
    #TAKING ONLY THE REAL PART FOR NOW!
    return Chi2DeltaDistXXZ_NR(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr)[0]*Lssp(npIR)

def Chi2DeltaDist_sps(omega,theta_center,omega_offset,capital_gamma,npVIS,Anr,Phi_nr):
    #TAKING ONLY THE REAL PART FOR NOW!
    return Chi2DeltaDistXZX_NR(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr)[0]*Lsps(npVIS)

def Chi2DeltaDist_pss(omega,theta_center,omega_offset,capital_gamma,npSF,Anr,Phi_nr):
    #TAKING ONLY THE REAL PART FOR NOW!
    return Chi2DeltaDistZXX_NR(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr)[0]*Lpss(npSF)


def CalculatedIntensity_PPP(omega,theta_center,omega_offset,capital_gamma,npSF,npVIS,npIR,ScaleFactor,Anr,Phi_nr):
    return ScaleFactor*abs(Chi2DeltaDist_ppp(omega,theta_center,omega_offset,capital_gamma,npSF,npVIS,npIR,Anr,Phi_nr))**2

def CalculatedIntensity_SSP(omega,theta_center,omega_offset,capital_gamma,npIR,ScaleFactor,Anr,Phi_nr):
    return ScaleFactor*abs(Chi2DeltaDist_ssp(omega,theta_center,omega_offset,capital_gamma,npIR,Anr,Phi_nr))**2

def CalculatedIntensity_PSS(omega,theta_center,omega_offset,capital_gamma,npVIS,ScaleFactor,Anr,Phi_nr):
    return ScaleFactor*abs(Chi2DeltaDist_pss(omega,theta_center,omega_offset,capital_gamma,npVIS,Anr,Phi_nr))**2

def CalculatedIntensity_SPS(omega,theta_center,omega_offset,capital_gamma,npSF,ScaleFactor,Anr,Phi_nr):
    return ScaleFactor*abs(Chi2DeltaDist_sps(omega,theta_center,omega_offset,capital_gamma,npSF,Anr,Phi_nr))**2


def XXZFomrulateDeltaDist(omega,theta_center,omega_offset,capital_gamma,ScaleFactorXXZ,Anr,Phi_nr):
    return abs(Chi2DeltaDistXXZ_NR_LabFrameTotal(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr,ScaleFactorXXZ))**2

def XZXFomrulateDeltaDist(omega,theta_center,omega_offset,capital_gamma,ScaleFactorXZX,Anr,Phi_nr):
    return abs(Chi2DeltaDistXZX_NR_LabFrameTotal(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr,ScaleFactorXZX))**2

def ZXXFomrulateDeltaDist(omega,theta_center,omega_offset,capital_gamma,ScaleFactorZXX,Anr,Phi_nr):
    return abs(Chi2DeltaDistZXX_NR_LabFrameTotal(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr,ScaleFactorZXX))**2

def ZZZFomrulateDeltaDist(omega,theta_center,omega_offset,capital_gamma,ScaleFactorZZZ,Anr,Phi_nr):
    return abs(Chi2DeltaDistZZZ_NR_LabFrameTotal(omega,theta_center,omega_offset,capital_gamma,Anr,Phi_nr,ScaleFactorZZZ))**2


OmegaVAL =  [i/2.0 for i in range(3100,3501)]
VSFGspec = [CalculatedIntensity_SSP(i,Theta_Center,Omega_offset,Capital_Gamma,\
                                    npIR,ScaleFactorXXZ,Anr,Phi_nr) for i in OmegaVAL]

VSFGspec2 = [CalculatedIntensity_SPS(i,Theta_Center,Omega_offset,Capital_Gamma,\
                                    npIR,ScaleFactorXXZ,Anr,Phi_nr) for i in OmegaVAL]
plt.figure()
plt.plot(OmegaVAL,VSFGspec)
plt.plot(OmegaVAL,VSFGspec2)
plt.xlabel('Wavelength')
plt.ylabel('Intensity (Normalized)')
plt.show()

