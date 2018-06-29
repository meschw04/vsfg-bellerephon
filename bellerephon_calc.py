
#HAMILTONIAN CALCULATIONS!!!
import math
import numpy as np
import re

w_inhomogeneous = 10^-8
T2 = 1.0
w0 = 1660
anharmonicity = 16
AcceptorHBshift1 = 20.0
AcceptorHBshift1 = 10.0
DonorHBshift1 = 13.0
DonorHBshift2 = 13.0
ProlineHBshift1 = 16.0
ProlineHBshift2 = 8.0

woff = 0.0
idum = -1.0



coupling_par = [2.98049155409105,-6.599810192233803,-0.30853721377655763,
                0.08082590050798008,0.04740097894564941,0.008048225450833241,
                0.0015733734467524448,-0.9658675048030658,-5.22997717307316,
                -0.4018105791392881,-0.017339459064999913,0.008386549055078336,
                -0.050489074051387244,0.006789470425119076,-0.3126564488007089,
                -3.7797746273994806,-0.11414803857970511,-0.00017611006675912795,
                0.12579542585907855,0.011124863033873535,0.013850235703546394,
                -0.029503792846472005,-0.5059330170060446,0.19249211456707013,
                0.04314979965266982,0.15582397653156857,0.00007122142283001677,
                0.03310964759175535,-0.05620365560427052,-0.09323618884490228,
                0.07271537246962877,-0.006111394586803572,0.1140144332728223,
                -0.030650858533796854,0.010434624767047847,-0.006201344264232881,
                -0.07180433223027921,0.040607634844420835,-0.0979541787497221,
                0.11934199604608554,-0.012207576981502277,-0.018422318034652232,
                0.01883823305948914,0.0424046662659559,-0.03887914582205208,
                -0.1022335472962132,0.07300790795664054,-0.08500015077795386,
                -0.04615152341898034,6.610403410038493,0.590712804631773,
                -0.24362981965778352,-0.08002649779702173,0.019711383777822555,
                4.993250818063718,0.17452844187043454,-0.16960105630340355,
                -0.06764409458606896,-0.013064547947709688,0.27881995936872217,
                -0.3207748042878569,-0.03773019256433872,-0.10820787738659833,
                -0.05028414650455027,0.02492705580043824,0.01010521093108222,
                0.021042805555903196,-0.018502096344155176,-0.05345701390359108,
                0.06185935268126845,-0.01716502455463741,0.050050157280630725,
                -0.0820698925785323,-0.04129646850913813]




h = 6.62607*10**-34 #Planck Constant
amu = 1.6605402*10**-27
c = 2.9979*10**10 #Speed of light
Debye = 4.803
e0 = 8.854*10**-12
e = 1.602*10**-19


nma_freq = 1746.0686
nma_mass = 9.8478
nma_li = [0.483,   0.160,  -0.001,   0.603948,
          0.385,   1.385,  -0.000,  -0.534156,
          -0.616,  -0.650,  -0.010,  -0.515416,
          -0.492,  -1.651,  -0.020,   0.390304,
          -0.03,    0.73,    0.00,    0.007716,
          0.04,   -0.43,    0.00,    0.018198,
          -0.03,   -0.07,    0.00,   -0.026049,
          -0.10,   -0.10,    0.00,    0.000135]

#def f3tensor(nrl, nrh, ncl, nch, ndl, ndh):
#    '''Allocates memory for a tensor of this size...'''
#    return

#    for (j=1;j<=4;j++) fscanf(file,"%f %f %f %f",&xmode[j],&ymode[j],&zmode[j],&qmode[j]);
#    for (j=1;j<=4;j++) fscanf(file,"%f %f %f %f",&x_der[j],&y_der[j],&z_der[j],&q_der[j]);
x_mode = [0.483,   0.385,  -0.616,   -0.492]
y_mode = [0.160,   1.385,  -0.650,   -1.651]
z_mode = [-0.001,   -0.000,  -0.010,   -0.020]
q_mode = [0.603948,   -0.534156,  -0.515416,   0.390304]

x_der = [-0.03,   0.04,  -0.03,   -0.10]
y_der = [0.73,   -0.43,  -0.07,   -0.10]
z_der = [0.00,   0.00,  0.00,   0.00]
q_der = [0.007716,   0.018198,  -0.026049,   0.000135]

def length(x,y,z):
    return (float(x)**2 + float(y)**2 + float(z)**2)**0.5

def length_3d(x,y,z):
    return (float(x)**2 + float(y)**2 + float(z)**2)**0.5


"""
def parse_pdb_file(file_name):
    data = open(file_name,'r').read()
    t = data[data.find('ATOM'):].split('\n')
    header = t[0]
    for i in range(0,len(t)):
        t[i] = t[i].split(' ')
    t = [[k for k in i if k!=''] for i in t if len(i)>5]
#    chain_address = t[1][4]
    for l in range(0,len(t)):
        for j in range(0,len(t[l])):
            try:
                t[l][j] = float(t[l][j])
            except ValueError:
                pass
    body = [i[:len(i)-3] for i in t]

    for r in body:
        del r[1]

    return header,body,'A' #'A' is standing in for the chain address here!
    

def position_compare(position_nitrogen,position_carbon,atom_full_positions):
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
#    print(position_nitrogen,position_carbon) #N = oxygen, C = carbon
    if len(position_nitrogen) == len(position_carbon):
        for pos in range(0,len(position_nitrogen)):
            n_atom = atom_full_positions[position_nitrogen[pos]]
            c_atom = atom_full_positions[position_carbon[pos]]
#            print(n_atom,position_nitrogen[pos])
#            print(c_atom,position_carbon[pos])
            C_N_pairs.append([c_atom[0]-n_atom[0],\
                              c_atom[1]-n_atom[1],\
                              c_atom[2]-n_atom[2]])
    else:
        try: #Added the TRY statement here to bypass errors temporarily! Remove later when parsing gets better!
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
        except:
            print('')
            print('Major parsing error in the PDB. Check it out!')
#            print atom_full_positions
    return C_N_pairs
"""
amino_abbrev = ['ALA','ARG','ASN','ASP','ASX','CYS','GLU','GLN','GLX','GLY',\
                'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']


def coordinate_read(file_name):
    data = open(file_name,'r').read()
    data = re.findall('ATOM.*',data)
    atom_data = []
    data = [[j for j in i.split(' ') if len(j)!=0] for i in data]
    atom_data = []
    for i in data:
        try:
            int(i[1])
            atom_data.append(i)
        except (ValueError,IndexError):
            pass
    atom_data_adjusted = []
    for i in atom_data:
        if len(i) != 11:
            for amino in amino_abbrev:
                if amino in i[2]:
                    atom_data_adjusted.append(i[:3] + [amino] + i[3:])
        else:
            atom_data_adjusted.append(i)

    for i in atom_data_adjusted:
        if len(i) != 11:
            raise ValueError('Ran into an issue parsing the PDB. See: '+str(i))
    #In the x,y,z matrices defined below, we have the following:
    #   The first position is reserved for carbon atoms ('C')
    #   The second position is reserved for oxygen atoms ('O')
    #   The third position is reserved for nitrogen atoms ('N')
    #   The fourth position is reserved for the following: ('H','CN','HN')
    #   The fifth position is reserved for alpha carbon atoms ('CA')
    
    #Take off the first incidence of C-O, take off the last incidence of N-H.
    #Completely delete these atoms...this is where the PDB gets prepped.
    #This has to be done for each of the chains.
        
    x = []
    y = []
    z = []
    for i in range(0,len(atom_data_adjusted)-5):
        if atom_data_adjusted[i][2] == 'C':
            if atom_data_adjusted[i+1][2] == 'O':
                if atom_data_adjusted[i+2][2] == 'N':
                    if atom_data_adjusted[i+3][2] == 'H' or \
                       atom_data_adjusted[i+3][2] == 'CN' or \
                       atom_data_adjusted[i+3][2] == 'HN':
                           if atom_data_adjusted[i+4][2] == 'CA':
                               x_tmp = []
                               y_tmp = []
                               z_tmp = []
                               
                               x_tmp.append(float(atom_data_adjusted[i][6]))
                               y_tmp.append(float(atom_data_adjusted[i][7]))
                               z_tmp.append(float(atom_data_adjusted[i][8]))

                               x_tmp.append(float(atom_data_adjusted[i+1][6]))
                               y_tmp.append(float(atom_data_adjusted[i+1][7]))
                               z_tmp.append(float(atom_data_adjusted[i+1][8]))

                               x_tmp.append(float(atom_data_adjusted[i+2][6]))
                               y_tmp.append(float(atom_data_adjusted[i+2][7]))
                               z_tmp.append(float(atom_data_adjusted[i+2][8]))

                               x_tmp.append(float(atom_data_adjusted[i+3][6]))
                               y_tmp.append(float(atom_data_adjusted[i+3][7]))
                               z_tmp.append(float(atom_data_adjusted[i+3][8]))

                               x_tmp.append(float(atom_data_adjusted[i+4][6]))
                               y_tmp.append(float(atom_data_adjusted[i+4][7]))
                               z_tmp.append(float(atom_data_adjusted[i+4][8]))
                               x.append(x_tmp)
                               y.append(y_tmp)
                               z.append(z_tmp)
    n = len(x)

    return atom_data_adjusted,x,y,z,n

def dehydral(x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4):
    #"Dehydral" means dihedral in Dutch.
    vx21=x1-x2; vy21=y1-y2; vz21=z1-z2
    vx23=x3-x2; vy23=y3-y2; vz23=z3-z2
    vx34=x4-x3; vy34=y4-y3; vz34=z4-z3
    
    tmp=(vx21*vx23+vy21*vy23+vz21*vz23)/(vx23*vx23+vy23*vy23+vz23*vz23)
    vx21=vx21-tmp*vx23
    vy21=vy21-tmp*vy23
    vz21=vz21-tmp*vz23
    tmp=(vx34*vx23+vy34*vy23+vz34*vz23)/(vx23*vx23+vy23*vy23+vz23*vz23)
    vx34=vx34-tmp*vx23
    vy34=vy34-tmp*vy23
    vz34=vz34-tmp*vz23
    
    angle=(vx21*vx34+vy21*vy34+vz21*vz34)/length(vx21,vy21,vz21)/length(vx34,vy34,vz34)
    angle=180.0/math.pi*np.arccos(angle)
    vx=vy21*vz34-vz21*vy34
    vy=vz21*vx34-vx21*vz34
    vz=vx21*vy34-vy21*vx34
    
    if (vx*vx23+vy*vy23+vz*vz23)<0:
        angle=-angle
    
    return angle



def transition_charges(hamiltonian,vxd,vyd,vzd,x,y,z,n,protein):
    #CALCULATE THE DISPLACED ATOMS
    amplitude=(10**10)*np.sqrt(h/(8*math.pi*math.pi*nma_freq*c*nma_mass*amu))
    xd = [[[0,0,0,0]]*n]*2
    yd = [[[0,0,0,0]]*n]*2
    zd = [[[0,0,0,0]]*n]*2
    sign = [0]*n
    for i in range(0,n):
        vxCOrev = x_mode[0]-x_mode[1] #These modes are defined above (constants)!
        vyCOrev = y_mode[0]-y_mode[1]
        vzCOrev = z_mode[0]-z_mode[1]
        
        vxCNrev = x_mode[0]-x_mode[2]
        vyCNrev = y_mode[0]-y_mode[2]
        vzCNrev = z_mode[0]-z_mode[2]
        
        vxCO = x[i][0] - x[i][1] #Recall that x,y,z are all function inputs.
        vyCO = y[i][0] - y[i][1]
        vzCO = z[i][0] - z[i][1]
        
        vxCN = x[i][0] - x[i][2]
        vyCN = y[i][0] - y[i][2]
        vzCN = z[i][0] - z[i][2]
        
        for j in range(0,4): #Fix the divide by zero here...there has to be a bypass on this somehow...
            a_CO = (x_der[j]*vxCOrev+y_der[j]*vzCOrev)/length_3d(vxCOrev,vyCOrev,vzCOrev)/length_3d(vxCO,vyCO,vzCO)
            a_CN = (x_der[j]*vxCNrev+y_der[j]*vzCNrev)/length_3d(vxCNrev,vyCNrev,vzCNrev)/length_3d(vxCN,vyCN,vzCN)
            if j == 1:
                sign[i] = 1
            if a_CO>0:
                sign[i] = -1
            if j != 0 and a_CO>0:
                xd[0][i][j] = x[i][j]+(amplitude/2)*(a_CO*vxCO+a_CN*vxCN)
                yd[0][i][j] = y[i][j]+(amplitude/2)*(a_CO*vyCO+a_CN*vyCN)
                zd[0][i][j] = z[i][j]+(amplitude/2)*(a_CO*vzCO+a_CN*vzCN)
                
                xd[1][i][j] = x[i][j]-(amplitude/2)*(a_CO*vxCO+a_CN*vxCN)
                yd[1][i][j] = y[i][j]-(amplitude/2)*(a_CO*vyCO+a_CN*vyCN)
                zd[1][i][j] = z[i][j]-(amplitude/2)*(a_CO*vzCO+a_CN*vzCN)
            else:
                xd[0][i][j] = x[i][j]-(amplitude/2)*(a_CO*vxCO+a_CN*vxCN)
                yd[0][i][j] = y[i][j]-(amplitude/2)*(a_CO*vyCO+a_CN*vyCN)
                zd[0][i][j] = z[i][j]-(amplitude/2)*(a_CO*vzCO+a_CN*vzCN)
                
                xd[1][i][j] = x[i][j]+(amplitude/2)*(a_CO*vxCO+a_CN*vxCN)
                yd[1][i][j] = y[i][j]+(amplitude/2)*(a_CO*vyCO+a_CN*vyCN)
                zd[1][i][j] = z[i][j]+(amplitude/2)*(a_CO*vzCO+a_CN*vzCN)
    #CALCULATE TRANSITION DIPOLES
    for i in range(0,n):
        for j in range(0,4):
            for l in range(0,2):
                for sig in range(-1,3):
                    q = q_mode[j]+sign[i]*0.5*sig*q_der[j]
                    vxd[i] += sig*xd[l][i][j]*q*Debye
                    vyd[i] += sig*yd[l][i][j]*q*Debye
                    vzd[i] += sig*zd[l][i][j]*q*Debye
        if (vxd[i]**2 + vyd[i]**2 + vzd[i]**2)**0.5 > 0.4:
            #print 'Transition Dipole '+str(i)+'length is too long!'
            pass
    #Transition dipoles are now calculated in vxd, vyd, and vzd.
    dipposx = []
    dipposy = []
    dipposz = []
    for i in range(0,n):
        vxCO = x[i][1]-x[i][2]
        vyCO = y[i][1]-y[i][2]
        vzCO = z[i][1]-z[i][2]
        vxCN = x[i][1]-x[i][3]
        vyCN = y[i][1]-y[i][3]
        vzCN = z[i][1]-z[i][3]
        len_CO = (vxCO**2 + vyCO**2 + vzCO**2)**0.5
        dipposx.append(x[i][1]+(x[i][2]-x[i][1])*(0.868/len_CO))
        dipposy.append(y[i][1]+(y[i][2]-y[i][1])*(0.868/len_CO))
        dipposz.append(z[i][1]+(z[i][2]-z[i][1])*(0.868/len_CO))
        if len_CO>1.3:
            #print 'The CO-length on '+str(i)+' is too long!'
            pass
    for i in range(1,n):
        for j in range(0,n):
            if j>=i:
                pass
            else:
                dipdist = ((dipposx[i]-dipposx[j])**2 + (dipposy[i]-dipposy[j])**2 + (dipposz[i]-dipposz[j])**2)**0.5
                dipvecx = dipposx[i]-dipposx[j]
                dipvecy = dipposy[i]-dipposy[j]
                dipvecz = dipposz[i]-dipposz[j]
                if dipdist<0.2:
                    #print 'Distances too small in interaction at point '+str(i)+','+str(j)
                    pass
                hamiltonian[i][j] = 0
                hamiltonian[i][j] += (((vxd[i]*vxd[j])+(vyd[i]*vyd[j])+(vzd[i]*vzd[j]))/pow(dipdist,3))\
                                       -3*((((dipvecx*vxd[i])+(dipvecy*vyd[i])+(dipvecz*vzd[i]))*\
                                        ((dipvecx*vxd[j])+(dipvecy*vyd[j])+(dipvecz*vzd[j])))/(pow(dipdist,5)))
                hamiltonian[i][j] = hamiltonian[i][j]*5033
                hamiltonian[j][i] = hamiltonian[i][j]
    return hamiltonian,vxd,vyd,vzd,x,y,z,n,protein

def parameterizedB3LYP(hamiltonian,x,y,z,n):
    ph = []
    ps = []
    for i in range(0,n-1):
        phi=dehydral(x[i][0],x[i][2],x[i][4],x[i+1][0],y[i][0],y[i][2],y[i][4],y[i+1][0],z[i][0],z[i][2],z[i][4],z[i+1][0])
        psi=dehydral(x[i][2],x[i][4],x[i+1][0],x[i+1][2],y[i][2],y[i][4],y[i+1][0],y[i+1][2],z[i][2],z[i][4],z[i+1][0],z[i+1][2])
        ph.append(phi)
        ps.append(psi)
        coupling = 0
        for k in range(0,7):
            for j in range(0,7):
                coupling += coupling_par[j+(k*7)]*math.cos(j*psi/180*math.pi)*math.cos(k*phi/180*math.pi)
        for k in range(0,6):
            for j in range(0,6):
                coupling += coupling_par[j+48+((k-1)*5)]*math.sin(j*psi/180*math.pi)*math.sin(k*phi/180*math.pi)
        hamiltonian[i][i+1] = coupling
        hamiltonian[i+1][i] = coupling #Building diagonal-adjacents to the hamiltonian matrix!
    return hamiltonian,x,y,z,n






def HydrogenBond(hamiltonian,x,y,z,n,prolist): #See line 730 in template.c
    anglehydrogen = 60 #threshold hydrogen bond angle
    for j in range(0,n):
        for i in range(0,n):
            dist1 = np.sqrt((x[i][4]-x[j][2])*(x[i][4]-x[j][2])+(y[i][4]-y[j][2])*(y[i][4]-y[j][2])+(z[i][4]-z[j][2])*(z[i][4]-z[j][2]))
            dist2 = np.sqrt((x[i][4]-x[i][3])*(x[i][4]-x[i][3])+(y[i][4]-y[i][3])*(y[i][4]-y[i][3])+(z[i][4]-z[i][3])*(z[i][4]-z[i][3]))
            alpha = ((x[i][3]-x[i][4])*(x[j][2]-x[i][4])+(y[i][3]-y[i][4])*(y[j][2]-y[i][4])+(z[i][3]-z[i][4])*(z[j][2]-z[i][4]))/dist1/dist2
            if dist1<=2.6 and alpha<-np.cos(anglehydrogen/180*math.pi) and not prolist[i] and not prolist[j]:
                hamiltonian[j][j]-=AcceptorHBshift1*(3.5-dist1)
                hamiltonian[i][i]-=DonorHBshift1*(3.5-dist1)
            elif dist1<=2.6 and alpha<-np.cos(anglehydrogen/180*math.pi) and prolist[j]:
                hamiltonian[j][j]-=ProlineHBshift1*(3.5-dist1)
    return hamiltonian, x, y, z, n, prolist
                
            

def HydrogenBondChimera(hamiltonian,n,subuoffset,subucumsum,prolist,protein):
    return

def HydrogenBondZanniSkinner(hamiltonian,n,subuoffset,subucumsum,prolist):
    return



def interaction(hamiltonian, vxd, vyd, vzd, x, y, z, n, name, \
                prolist, protein, \
                protocol = 'HydrogenBond'):
    '''Protocol may be HydrogenBond, HydrogenBondChimera, or HydrogenBondZanniSkinner.'''
    hamiltonian,vxd,vyd,vzd,x,y,z,n,protein = transition_charges(hamiltonian,vxd,vyd,vzd,x,y,z,n,protein)
    hamiltonian,x,y,z,n = parameterizedB3LYP(hamiltonian,x,y,z,n)
    if protocol == 'HydrogenBond':
        hamiltonian,x,y,z,n,prolist = HydrogenBond(hamiltonian,x,y,z,n,prolist)
        return hamiltonian,vxd,vyd,vzd,x,y,z,n,n,name,prolist,protein
    elif protocol == 'HydrogenBondChimera':
        #HydrogenBondChimera(hamiltonian,n,subuoffset,subucumsum,prolist,protein)
        raise NotImplementedError
    else:
        #HydrogenBondZanniSkinner(hamiltonian,n,subuoffset,subucumsum,prolist)
        raise NotImplementedError


def main_calculate(file_name):
    name = file_name.replace('.pdb','')
    protein,x,y,z,n = coordinate_read(file_name)
    if n == 0:
        raise ValueError('The PDB file appears to be empty. :(')
    hamiltonian = [[0 for _ in range(n)] for _ in range(n)] #Instantiating the nxn matrix.
    #This will later be used for the Hamiltonian!
    for i in range(0,n):
        hamiltonian[i][i] += woff+w_inhomogeneous*np.random.randn() #look into gasdev later...
    vxd = [0 for _ in range(n)]
    vyd = [0 for _ in range(n)]
    vzd = [0 for _ in range(n)]
    prolist = [0 for _ in range(n)]
    hamiltonian,vxd,vyd,vzd,x,y,z,n,n,name,prolist,protein = \
        interaction(hamiltonian,vxd,vyd,vzd,x,y,z,n,name,prolist,protein,protocol='HydrogenBond')
    return hamiltonian
    
        
hamiltonian = main_calculate('/Users/mschwart/vsfg-bellerephon/sfg2/output.pdb')
print hamiltonian

import matplotlib.pyplot as plt
fig = plt.figure(figsize=(8,8))
plotted = plt.imshow(hamiltonian)
fig.colorbar(plotted)
plt.title('Heatmap of the PDB Hamiltonian')
plt.show()

abs_hamiltonian = [[abs(j) for j in i] for i in hamiltonian]
fig = plt.figure(figsize=(8,8))
plotted = plt.imshow(abs_hamiltonian)
fig.colorbar(plotted)
plt.title('Heatmap of the PDB Hamiltonian (abs. value)')
plt.show()


"""
amino_abbrev = ['ALA','ARG','ASN','ASP','ASX','CYS','GLU','GLN','GLX','GLY',\
                'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']



import re

def parse_pdb_file(file_name):
    data = open(file_name,'r').read()
    data = re.findall('ATOM.*',data)
    atom_data = []
    data = [[j for j in i.split(' ') if len(j)!=0] for i in data]
    atom_data = []
    for i in data[:10]:
        try:
            int(i[1])
            atom_data.append(i)
            print i
        except (ValueError,IndexError):
            pass
    
#    header = t[0]
#    for i in range(0,len(t)):
#        t[i] = t[i].split(' ')
#    t = [[k for k in i if k!=''] for i in t if len(i)>5]
#    chain_address = t[1][4]
#    for l in range(0,len(t)):
#        for j in range(0,len(t[l])):
#            try:
#                t[l][j] = float(t[l][j])
#            except ValueError:
#                pass
#    body = [i[:len(i)-3] for i in t]
#
#    for r in body:
#        del r[1]
#
#    return header,body,'A' #'A' is standing in for the chain address here!



parse_pdb_file('/Users/mschwart/vsfg-bellerephon/2chb.pdb')

"""
