
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

x_mode = [0.483,   0.385,  -0.616,   -0.492]
y_mode = [0.160,   1.385,  -0.650,   -1.651]
z_mode = [-0.001,   -0.000,  -0.010,   -0.020]
q_mode = [0.603948,   -0.534156,  -0.515416,   0.390304]

x_der = [-0.03,   0.04,  -0.03,   -0.10]
y_der = [0.73,   -0.43,  -0.07,   -0.10]
z_der = [0.00,   0.00,  0.00,   0.00]
q_der = [0.007716,   0.018198,  -0.026049,   0.000135]

amino_abbrev = ['ALA','ARG','ASN','ASP','ASX','CYS','GLU','GLN','GLX','GLY',\
                'HIS','ILE','LEU','LYS','MET','PHE','PRO','SER','THR','TRP','TYR','VAL']

def length(x,y,z):
    return (float(x)**2 + float(y)**2 + float(z)**2)**0.5

def length_3d(x,y,z):
    return (float(x)**2 + float(y)**2 + float(z)**2)**0.5



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
    xd1 = [[0 for _ in range(0,4)] for _ in range(0,n)]
    yd1 = [[0 for _ in range(0,4)] for _ in range(0,n)]
    zd1 = [[0 for _ in range(0,4)] for _ in range(0,n)]
    xd2 = [[0 for _ in range(0,4)] for _ in range(0,n)]
    yd2 = [[0 for _ in range(0,4)] for _ in range(0,n)]
    zd2 = [[0 for _ in range(0,4)] for _ in range(0,n)]
    xd = [xd1,xd2]
    yd = [yd1,yd2]
    zd = [zd1,zd2]
    #vzd = [0 for _ in range(n)]    
    #    yd = [[[0,0,0,0]]*n]*2
    #zd = [[[0,0,0,0]]*n]*2
    sign = [0 for _ in range(n)]
    dipposx = []
    dipposy = []
    dipposz = []
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
        
        for j in range(0,4): #x_der is a constant, as as y_der
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
#                xd[0][i][j] = x[i][j]+(amplitude/2)*(a_CO*vxCO+a_CN*vxCN)
#                yd[0][i][j] = y[i][j]+(amplitude/2)*(a_CO*vyCO+a_CN*vyCN)
#                zd[0][i][j] = z[i][j]+(amplitude/2)*(a_CO*vzCO+a_CN*vzCN)
#                xd[1][i][j] = x[i][j]+(amplitude/2)*(a_CO*vxCO+a_CN*vxCN)
#                yd[1][i][j] = y[i][j]+(amplitude/2)*(a_CO*vyCO+a_CN*vyCN)
#                zd[1][i][j] = z[i][j]+(amplitude/2)*(a_CO*vzCO+a_CN*vzCN)
#                if (amplitude/2)*(a_CO*vzCO+a_CN*vzCN) > 0.0001:
#                    print (amplitude/2)*(a_CO*vzCO+a_CN*vzCN)
#                    print amplitude
            else:
#                xd[0][i][j] = x[i][j]-(amplitude/2)*(a_CO*vxCO+a_CN*vxCN)
#                yd[0][i][j] = y[i][j]-(amplitude/2)*(a_CO*vyCO+a_CN*vyCN)
#                zd[0][i][j] = z[i][j]-(amplitude/2)*(a_CO*vzCO+a_CN*vzCN)                
#                xd[1][i][j] = x[i][j]-(amplitude/2)*(a_CO*vxCO+a_CN*vxCN)
#                yd[1][i][j] = y[i][j]-(amplitude/2)*(a_CO*vyCO+a_CN*vyCN)
#                zd[1][i][j] = z[i][j]-(amplitude/2)*(a_CO*vzCO+a_CN*vzCN)
                xd[0][i][j] = x[i][j]-(amplitude/2)*(a_CO*vxCO+a_CN*vxCN)
                yd[0][i][j] = y[i][j]-(amplitude/2)*(a_CO*vyCO+a_CN*vyCN)
                zd[0][i][j] = z[i][j]-(amplitude/2)*(a_CO*vzCO+a_CN*vzCN)                
                xd[1][i][j] = x[i][j]+(amplitude/2)*(a_CO*vxCO+a_CN*vxCN)
                yd[1][i][j] = y[i][j]+(amplitude/2)*(a_CO*vyCO+a_CN*vyCN)
                zd[1][i][j] = z[i][j]+(amplitude/2)*(a_CO*vzCO+a_CN*vzCN)
#    CALCULATE TRANSITION DIPOLES
            for sig in range(-1,3):
                q = q_mode[j]+sign[i]*0.5*sig*q_der[j]
#                vxd[i] = sig*(xd[0][i][j]+xd[1][i][j])*q*Debye
#                vyd[i] = sig*(yd[0][i][j]+yd[1][i][j])*q*Debye
#                vzd[i] = sig*(zd[0][i][j]+zd[1][i][j])*q*Debye
                vxd[i] = sig*(xd[0][i][j]-xd[1][i][j])*q*Debye #Roeters used plus signs here.
                vyd[i] = sig*(yd[0][i][j]-yd[1][i][j])*q*Debye #We subtracted instead by ignoring
                vzd[i] = sig*(zd[0][i][j]-zd[1][i][j])*q*Debye #the "sign" list instantiated above.
#                if vxd[i] != 0 and vyd[i] !=0 and vzd[i] != 0:
#                    print vxd[i]
#                    print vyd[i]
#                    print vzd[i]
#                print (vxd[i]**2 + vyd[i]**2 + vzd[i]**2)**0.5
#        print [vxd[i],vyd[i],vzd[i]]
#        print (vxd[i]**2 + vyd[i]**2 + vzd[i]**2)**0.5
        if (vxd[i]**2 + vyd[i]**2 + vzd[i]**2)**0.5 > 0.4:
            #print 'Transition Dipole '+str(i)+'length is too long!'
            pass
        #Transition dipoles are now calculated in vxd, vyd, and vzd.
        vxCO = x[i][1]-x[i][2]
        vyCO = y[i][1]-y[i][2]
        vzCO = z[i][1]-z[i][2]
        vxCN = x[i][1]-x[i][3]
        vyCN = y[i][1]-y[i][3]
        vzCN = z[i][1]-z[i][3]
        len_CO = (vxCO**2 + vyCO**2 + vzCO**2)**0.5
#        print 'len_CO = '+str(len_CO)
        dipposx.append(x[i][1]+(x[i][2]-x[i][1])*(0.868/len_CO))
        dipposy.append(y[i][1]+(y[i][2]-y[i][1])*(0.868/len_CO))
        dipposz.append(z[i][1]+(z[i][2]-z[i][1])*(0.868/len_CO))
        if len_CO>1.3:
            #print 'The CO-length on '+str(i)+' is too long!'
            pass
#    print x[:10]
#    print y[:10]
#    print z[:10]
#    print sign
#    print xd[:2]
#    print yd[:2]
#    print zd[:2]
    dipdist2 = []
    for i in range(1,n):#Note: This is in range 1 to n, not 0 to n!
        dipdist1 = []
        for j in range(0,n):
            if j>=i:
                pass
            else:
                dipdist = ((dipposx[i]-dipposx[j])**2 + (dipposy[i]-dipposy[j])**2 + (dipposz[i]-dipposz[j])**2)**0.5
                dipdist1.append(dipdist)
                dipvecx = dipposx[i]-dipposx[j]
                dipvecy = dipposy[i]-dipposy[j]
                dipvecz = dipposz[i]-dipposz[j]
                if dipdist<0.2:
                    #print 'Distances too small in interaction at point '+str(i)+','+str(j)
                    pass
                hamiltonian[i][j] = 0
                hamiltonian[i][j] += (((vxd[i]*vxd[j])+(vyd[i]*vyd[j])+(vzd[i]*vzd[j]))/(dipdist**3))\
                                       -3*((((dipvecx*vxd[i])+(dipvecy*vyd[i])+(dipvecz*vzd[i]))*\
                                        ((dipvecx*vxd[j])+(dipvecy*vyd[j])+(dipvecz*vzd[j])))/((dipdist**5)))
                hamiltonian[i][j] = hamiltonian[i][j]*5033
                hamiltonian[j][i] = hamiltonian[i][j] #Where the non-coupled diagonals are calculated...
        dipdist2.append(dipdist1)
    print dipdist2
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
            for j in range(0,6): #Check the nearest neighbor couplings...
                coupling += coupling_par[j+48+((k-1)*5)]*math.sin(j*psi/180*math.pi)*math.sin(k*phi/180*math.pi)
        hamiltonian[i][i+1] = coupling #This is the nearest neighbor coupling
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
    for i in range(0,n):
        for j in range(0,n):
            if i == j:
                hamiltonian[i][j] = 0

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
#print hamiltonian

import matplotlib.pyplot as plt
fig = plt.figure(figsize=(4,4))

l = []
for k in hamiltonian[:25]:
    l.append(k[:25])

plotted = plt.imshow(l)
fig.colorbar(plotted)
plt.title('Heatmap of the PDB Hamiltonian')
plt.show()

#abs_hamiltonian = [[abs(j) for j in i] for i in hamiltonian]
#fig = plt.figure(figsize=(8,8))
#plotted = plt.imshow(abs_hamiltonian)
#fig.colorbar(plotted)
#plt.title('Heatmap of the PDB Hamiltonian (abs. value)')
#plt.show()


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
qrs = l

haa_file = '/Users/mschwart/vsfg-bellerephon/sfg2/output_prepped_haa.txt'
t = open(haa_file,'rb').read()
q = t.split('\n')
p = []
for i in q:
    j = []
    for l in i.split(' '):
        try:
            j.append(float(l))
        except:
            pass
    p.append(j)
p = p[:len(p)-1]

o = []
for i in p[:25]:
    o.append(i[:25])

fig = plt.figure(figsize=(4,4))
plotted = plt.imshow(o)
fig.colorbar(plotted)
plt.title('Heatmap of the PDB Hamiltonian (abs. value)')
plt.show()

#print o
#print l

#print o
#print ''
#print qrs

#p1 = o[0][2:25]
#print ''
#p2 = qrs[0][2:25]

#print [i/j for (i,j) in zip(p1,p2)]


'''
7/19/2018
Vance and I worked on the stuff. We got the following outputs from Roeters code.
The first block represents the dipole distances, which are different from ours for some reason.
The second block represents the outputs vxd, vyd, vzd, etc., see line 624 of Roeters.
The main issue here is that the k-nearest neighbors (not the diagonals or diagonal adjacents)
are having major issues. This is what needs to be debugged; that is, all the Hamiltonian
values that aren't a diagonal or adjacent to a diagonal. The next block is another look at
the various dipole distances in the system. We just gotta fix this issue then we'll be groovy.

first group of the dipole distances...
nan
  nan
3.050
  nan
4.926
3.518
  nan
4.857
4.901
3.628
  nan
6.086
4.949
5.071
3.353
  nan
8.519
6.507
5.044
5.135
3.720
  nan
8.418
7.533
5.600
3.747
3.838
3.299
  nan
11.955
10.887
9.201
7.392
6.435
5.233
3.780
  nan
12.762
10.873
9.452
8.851
7.046
4.496
5.511
3.620
  nan
16.045
14.009
12.540
12.152
10.327
7.647
8.696
6.066
3.329
  nan
18.056
16.399
14.209
13.586
12.653
9.922
9.873
6.898
6.077
4.226
  nan
18.993
17.718
16.017
14.466
13.207
11.520
10.812
7.079
7.424
5.971
4.370
  nan
23.143
21.827
20.135
18.624
17.306
15.571
14.971
11.240
11.298
9.124
7.098
4.163
  nan
297.513
295.511
293.364
293.084
291.872
289.089
289.356
286.012
284.897
281.685
279.510
279.256
275.340
  nan
296.708
294.692
292.578
292.296
291.048
288.275
288.566
285.210
284.072
280.851
278.718
278.445
274.523
3.490
  nan
299.813
297.814
295.710
295.389
294.133
291.388
291.655
288.281
287.169
283.952
281.815
281.492
277.553
5.711
4.513
  nan
297.745
295.767
293.673
293.307
292.044
289.331
289.569
286.175
285.094
281.882
279.742
279.361
275.404
8.413
7.084
4.295
  nan
298.949
296.995
294.914
294.497
293.226
290.550
290.757
287.340
286.293
283.087
280.944
280.499
276.523
12.463
11.405
8.001
4.339
  nan
302.795
300.846
298.749
298.336
297.082
294.398
294.597
291.186
290.148
286.945
284.783
284.350
280.378
12.701
12.379
8.147
6.389
4.193
  nan
305.431
303.495
301.399
300.963
299.711
297.042
297.223
293.803
292.785
289.587
287.416
286.957
282.978
15.489
15.398
11.039
9.638
6.937
3.253





x,y,z & length of dipole 1:    nan    nan    nan    nan
Distance between C-atom of backbone and the dipole at peptide bond 1: nan
Transition Dipole peptide bond 1: nan nan nan, Length: nan and CO vector and length: (0.000000 0.000000 0.000000, 0.000000)

 x,y,z & length of dipole 2:  0.196 -0.187 -0.240  0.362
Distance between C-atom of backbone and the dipole at peptide bond 2: 0.868000
Transition Dipole peptide bond 2: 0.195559 -0.187031 -0.240357, Length: 0.361932 and CO vector and length: (-0.598007 1.041992 0.259998, 1.229211)

 x,y,z & length of dipole 3:  0.335 -0.112  0.112  0.371
Distance between C-atom of backbone and the dipole at peptide bond 3: 0.868000
Transition Dipole peptide bond 3: 0.335226 -0.111812 0.112163, Length: 0.370754 and CO vector and length: (-0.871994 0.864014 -0.038002, 1.228144)

 x,y,z & length of dipole 4: -0.064 -0.339  0.127  0.367
Distance between C-atom of backbone and the dipole at peptide bond 4: 0.868000
Transition Dipole peptide bond 4: -0.063739 -0.338550 0.127362, Length: 0.367287 and CO vector and length: (-0.222000 1.207977 0.067997, 1.230088)

 x,y,z & length of dipole 5:  0.038 -0.255 -0.281  0.381
Distance between C-atom of backbone and the dipole at peptide bond 5: 0.868000
Transition Dipole peptide bond 5: 0.037926 -0.254789 -0.281348, Length: 0.381461 and CO vector and length: (-0.622009 0.854004 0.625999, 1.228045)

 x,y,z & length of dipole 6:  0.342 -0.057 -0.134  0.371
Distance between C-atom of backbone and the dipole at peptide bond 6: 0.868000
Transition Dipole peptide bond 6: 0.341585 -0.056926 -0.133761, Length: 0.371232 and CO vector and length: (-0.972000 0.735992 0.155998, 1.229148)

 x,y,z & length of dipole 7:  0.131 -0.261  0.227  0.370
Distance between C-atom of backbone and the dipole at peptide bond 7: 0.868000
Transition Dipole peptide bond 7: 0.130840 -0.260933 0.227378, Length: 0.370008 and CO vector and length: (-0.209991 1.186005 -0.243999, 1.228918)



Dipoles:  2   1 Difference vector:   nan   nan   nan
Difference vector length calc. 1:   nan
Dipole length:  2 = 449.049,  1 =   nan
Difference vector length calc. 2:   nan
Dipoles:  3   1 Difference vector:   nan   nan   nan
Difference vector length calc. 1:   nan
Dipole length:  3 = 446.620,  1 =   nan
Difference vector length calc. 2:   nan
Dipoles:  3   2 Difference vector: -0.656 -2.153 2.058
Difference vector length calc. 1: 3.050
Dipole length:  3 = 446.620,  2 = 449.049
Difference vector length calc. 2: 3.050
Dipoles:  4   1 Difference vector:   nan   nan   nan
Difference vector length calc. 1:   nan
Dipole length:  4 = 444.188,  1 =   nan
Difference vector length calc. 2:   nan
Dipoles:  4   2 Difference vector: -1.880 -4.549 -0.208
Difference vector length calc. 1: 4.926
Dipole length:  4 = 444.188,  2 = 449.049
Difference vector length calc. 2: 4.926
Dipoles:  4   3 Difference vector: -1.224 -2.396 -2.267
Difference vector length calc. 1: 3.518
Dipole length:  4 = 444.188,  3 = 446.620
Difference vector length calc. 2: 3.518
Dipoles:  5   1 Difference vector:   nan   nan   nan
Difference vector length calc. 1:   nan
Dipole length:  5 = 446.020,  1 =   nan
Difference vector length calc. 2:   nan
Dipoles:  5   2 Difference vector: 1.321 -4.274 -1.893
Difference vector length calc. 1: 4.857
Dipole length:  5 = 446.020,  2 = 449.049
Difference vector length calc. 2: 4.857
Dipoles:  5   3 Difference vector: 1.977 -2.121 -3.951
Difference vector length calc. 1: 4.901
Dipole length:  5 = 446.020,  3 = 446.620
Difference vector length calc. 2: 4.901
Dipoles:  5   4 Difference vector: 3.201 0.275 -1.684
Difference vector length calc. 1: 3.628
Dipole length:  5 = 446.020,  4 = 444.188
'''

