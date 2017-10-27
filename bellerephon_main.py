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
    C_N_pairs = []
    if len(position_nitrogen) == len(position_carbon):
        for pos in range(0,len(position_nitrogen)):
            n_atom = atom_full_positions[position_nitrogen[pos]]
            c_atom = atom_full_positions[position_carbon[pos]]
            C_N_pairs.append([c_atom[0]-n_atom[0],\
                              c_atom[1]-n_atom[1],\
                              c_atom[2]-n_atom[2]])
    else:
        for i in range(0,len(max(position_nitrogen,position_carbon))):
#            print(position_nitrogen)
#            print(position_carbon)
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


header, body, chain_address = parse_pdb_file('C:\Users\Marcus\Documents\Fall 2017\sfg2\output_CHAIN_A.pdb')

#header, body, chain_address = parse_pdb_file('C:\Users\Marcus\Documents\Fall 2017\sfg2\output.pdb')

#header, body, chain_address = parse_pdb_file('C:\Users\Marcus\Documents\Fall 2017\sfg2\output_TZ.pdb')


#This is the cR from the Mathematica code
trans_dip_mom = pd.read_csv('C:\Users\Marcus\Documents\Fall 2017\sfg2\output_prepped_transdipmom_Hamm.txt',header=None)
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



eigenvalues = pd.read_csv('C:\Users\Marcus\Documents\Fall 2017\sfg2\output_prepped_eval.txt',header=None)
eig_val = []
for i in eigenvalues[0].tolist():
    eig_val.append(float(i)+1660)
eigenvalues = eig_val
eigenvectors = pd.read_csv('C:\Users\Marcus\Documents\Fall 2017\sfg2\output_prepped_evec.txt',header=None)

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
    C_N_pairs = position_compare(position_nitrogen,position_carbon)
    C_O_pairs = position_compare(position_oxygen,position_carbon)

    c = []
    for i in C_O_pairs:
        norm_int = (i[0]**2 + i[1]**2 + i[2]**2)**0.5
        new_norm = []
        for k in i:
            new_norm.append(k/norm_int)
        c.append(new_norm)
    
    a = []
    for i in range(0,len(c)):
        numer = list(np.cross(c[i],C_N_pairs[i]))
        denom = (numer[0]**2+numer[1]**2+numer[2]**2)**0.5
        inter = []
        for j in numer:
            inter.append(j/denom)
        a.append(inter)
    
    
    b = []
    for i in range(0,len(c)):
        crossed = list(np.cross(c[i],a[i]))
        b.append(crossed)
    
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
    
    #Iir1mode=Table[ Sum[EigenVectorsHam[[i,j]]cR[[j]],{j,LH}],{i,LH}];
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
    
    
    #Iir1modeAbs=Table[(Norm[Sum[ EigenVectorsHam[[i,j]]cR[[j]],{j,LH}]])^2,{i,LH}];
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

#IRspec=Flatten[Table[IRspectrum[\[Omega]],{\[Omega],1550,1750,0.5}]];
#OmegaVal=Table[\[Omega],{\[Omega],1550,1750,0.5}];
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
### LATER ON, PIPE OUT A CSV WITH OMEGAVAL, IRSPEC

#Iraman1mode=Table[Sum[EigenVectorsHam[[i,j]]RamanTensor[[j]],{j, LH}],{i,LH}];
#Iraman1modeIsotropic=Table[1/9 (Tr[Sum[EigenVectorsHam[[i,j]]RamanTensor[[j]],{j, LH}]])^2 ,{i,LH}];


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

#Iraman1modeAnisotropic=Table[ 1/2 (  ( Iraman1mode[[i,1,1]]- Iraman1mode[[i,2,2]] )^2  + (Iraman1mode[[i,2,2]]- Iraman1mode[[i,3,3]])^2 +( Iraman1mode[[i,3,3]]- Iraman1mode[[i,1,1]] )^2)+
#3/4 (  ( Iraman1mode[[i,1,2]]- Iraman1mode[[i,2,1]] )^2  + (Iraman1mode[[i,2,3]]- Iraman1mode[[i,3,2]])^2 +( Iraman1mode[[i,3,1]]- Iraman1mode[[i,1,3]] )^2 ),{i,LH}];
#
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

#RamanspectrumIsotropic[\[Omega]_,\[CapitalGamma]Raman_]:=Sum[(Iraman1modeIsotropic[[i]](\[CapitalGamma]Raman/Pi))/((\[Omega]-EigenValuesHam[[i]]-\[Omega]offsetIR)^2+ \[CapitalGamma]Raman^2),{i,LH}]
#RamanspectrumAnisotropic[\[Omega]_,\[CapitalGamma]Raman_]:=Sum[(Iraman1modeAnisotropic[[i]](\[CapitalGamma]Raman/Pi))/((\[Omega]-EigenValuesHam[[i]]-\[Omega]offsetIR)^2+ \[CapitalGamma]Raman^2),{i,LH}]

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
    

Capital_Gamma_Raman = 9.0

#Isfg1mode=Table[Iraman1mode[[q,l,m]]Iir1mode[[q,n]],{q,LH},{l,3},{m,3},{n,3}];
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

    
    
    returned_val = [Ns*i for i in returned_val]
    return returned_val




Chi2DeltaDistZZZEnsembleDeltaDist(4)