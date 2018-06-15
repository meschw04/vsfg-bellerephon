
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
'''
/* Need this one! */
void transition_charges(float **H1,float *vxd,float *vyd,float *vzd,float **x,float **y,float **z,int n, int istruct, const char *protein)
{
    char temp[100]; #DONE
    float xmode[5],ymode[5],zmode[5],qmode[5],Freq,Mass; #DONE
    float x_der[5],y_der[5],z_der[5],q_der[5]; #DONE
    float vxCOrev,vyCOrev,vzCOrev,vxCO,vyCO,vzCO; #DONE
    float *dipposx,*dipposy,*dipposz; #DONE
    float vxCNrev,vyCNrev,vzCNrev,vxCN,vyCN,vzCN; #DONE
    float aCO,aCN,Amplitude,xtmp,q,sig,q1,q2,q3,q4,r1,r2,r3,r4; #DONE
    float ***xd, ***yd, ***zd; #DONE
    FILE *file,*transdipmomfile; #DONE
    int i,i1,i2,j,j1,j2,l,sign[nmax]; #DONE
    
    xd=f3tensor(1,2,1,n,1,4);yd=f3tensor(1,2,1,n,1,4);zd=f3tensor(1,2,1,n,1,4); #DONE
    dipposx=vector(1,n);dipposy=vector(1,n);dipposz=vector(1,n); #DONE

    
    /*Load parameter files for displacements, charges and charge flows*/
    if((file=fopen("nma.par","rt"))==NULL) {printf ("\nCan't open %s\n","nma.par");exit(1);}
    fscanf(file,"%f ",&Freq);
    fscanf(file,"%f ",&Mass);
    for (j=1;j<=4;j++) fscanf(file,"%f %f %f %f",&xmode[j],&ymode[j],&zmode[j],&qmode[j]);
    for (j=1;j<=4;j++) fscanf(file,"%f %f %f %f",&x_der[j],&y_der[j],&z_der[j],&q_der[j]);
    
    
    /*Calculate displaced atoms*/ ###
    Amplitude=1e10*sqrt(h/(8*PI*PI*Freq*c*Mass*amu)); ###
    for (i=1;i<=n;i++)
    {
        vxCOrev=xmode[1]-xmode[2];vyCOrev=ymode[1]-ymode[2];vzCOrev=zmode[1]-zmode[2];
        vxCNrev=xmode[1]-xmode[3];vyCNrev=ymode[1]-ymode[3];vzCNrev=zmode[1]-zmode[3];
        vxCO=x[i][1]-x[i][2];vyCO=y[i][1]-y[i][2];vzCO=z[i][1]-z[i][2];
        vxCN=x[i][1]-x[i][3];vyCN=y[i][1]-y[i][3];vzCN=z[i][1]-z[i][3];
        
        for (j=1;j<=4;j++)
        {
            aCO=(x_der[j]*vxCOrev+y_der[j]*vyCOrev+z_der[j]*vzCOrev)/length(vxCOrev,vyCOrev,vzCOrev)/length(vxCO,vyCO,vzCO);
            aCN=(x_der[j]*vxCNrev+y_der[j]*vyCNrev+z_der[j]*vzCNrev)/length(vxCNrev,vyCNrev,vzCNrev)/length(vxCN,vyCN,vzCN);
            if (j==1)
            {
                sign[i]=1;                                /*Test sign of Eigenmode by movement of C-atom)*/
                if (aCO>0) sign[i]=-1;
            }
            
            xd[1][i][j]=x[i][j]-sign[i]*Amplitude/2.*(aCO*vxCO+aCN*vxCN);
            yd[1][i][j]=y[i][j]-sign[i]*Amplitude/2.*(aCO*vyCO+aCN*vyCN);
            zd[1][i][j]=z[i][j]-sign[i]*Amplitude/2.*(aCO*vzCO+aCN*vzCN);
            xd[2][i][j]=x[i][j]+sign[i]*Amplitude/2.*(aCO*vxCO+aCN*vxCN);
            yd[2][i][j]=y[i][j]+sign[i]*Amplitude/2.*(aCO*vyCO+aCN*vyCN);
            zd[2][i][j]=z[i][j]+sign[i]*Amplitude/2.*(aCO*vzCO+aCN*vzCN);
        }
    }
    
  
  /*First calculate transition dipoles*/
    for (i=1;i<=n;i++)
    {
        vxd[i]=0;vyd[i]=0;vzd[i]=0;
        for (j=1;j<=4;j++)
            for (l=1,sig=-1;l<=2;l++,sig+=2)
            {
                q=qmode[j]+sign[i]*0.5*sig*q_der[j];
                vxd[i]+=sig*xd[l][i][j]*q*Debye;
                vyd[i]+=sig*yd[l][i][j]*q*Debye;
                vzd[i]+=sig*zd[l][i][j]*q*Debye;
            }
        printf("Transition Dipoles: %f %f %f\n",vxd[i],vyd[i],vzd[i]);
        printf("Transition Dipoles Length: %f\n",length(vxd[i],vyd[i],vzd[i]));
       if (length(vxd[i],vyd[i],vzd[i])>0.4) {
            printf("Transition dipole (i is %d) length is too long! Coordinates are: C - %f %f %f, N - %f %f %f and O-  %f %f %f", i, x[i][1], y[i][1], z[i][1], x[i][3], y[i][3], z[i][3], x[i][2], y[i][2], z[i][2]);        }

    
    }     
    
    sprintf(temp,"%s_transdipmom_Hamm.txt",protein);
    transdipmomfile = fopen(temp,"wt");
   
    /* Print transition dipole moments to file: */
    for(i=1;i<=n;i++){
        fprintf(transdipmomfile," %7.2lf %7.2lf %7.2lf \n",vxd[i],vyd[i],vzd[i]);
    }
    
            fclose(transdipmomfile);
   
    /*Calculate transition dipole positions (by Steven Roeters)*/
    for (i=1;i<=n;i++)
    {
        // Construction of CO and CN-vectors:
        vxCO=x[i][1]-x[i][2];vyCO=y[i][1]-y[i][2];vzCO=z[i][1]-z[i][2];
        vxCN=x[i][1]-x[i][3];vyCN=y[i][1]-y[i][3];vzCN=z[i][1]-z[i][3];
        
        // Placement of these vectors in space, according to Torii et al. '91, Krimm et al. '72 and Sandeman '55:
        dipposx[i] = x[i][1] + ( (x[i][2]-x[i][1]) * ( 0.868 / length(vxCO,vyCO,vzCO) ) );
        dipposy[i] = y[i][1] + ( (y[i][2]-y[i][1]) * ( 0.868 / length(vxCO,vyCO,vzCO) ) );
        dipposz[i] = z[i][1] + ( (z[i][2]-z[i][1]) * ( 0.868 / length(vxCO,vyCO,vzCO) ) );
                  printf("\n x,y,z & length of dipole %d: %6.3f %6.3f %6.3f %6.3f\n",i,vxd[i],vyd[i],vzd[i], length(vxd[i],vyd[i],vzd[i]));
    
            printf("Distance between C-atom of backbone and the dipole at peptide bond %d: %f\n",i,length((x[i][2]-x[i][1])*(0.868/length(vxCO,vyCO,vzCO)),(y[i][2]-y[i][1])*(0.868/length(vxCO,vyCO,vzCO)),(z[i][2]-z[i][1])*(0.868/length(vxCO,vyCO,vzCO))));
        
 printf("Transition Dipole peptide bond %d: %f %f %f, Length: %f and CO vector and length: (%f %f %f, %f)\n",i,vxd[i],vyd[i],vzd[i],length(vxd[i],vyd[i],vzd[i]),vxCO,vyCO,vzCO,length(vxCO,vyCO,vzCO));
        
        if (length(vxCO,vyCO,vzCO)>1.3) {
                    printf("The CO-length is too long (%f A)!\n",length(vxCO,vyCO,vzCO));
        }

    }
      
    {
    int i;
    float dipdist,dipvecx,dipvecy,dipvecz;
for (i=2;i<=n;i++)
for (j=1;j<i;j++)
{
            dipdist=sqrt((dipposx[i]-dipposx[j])*(dipposx[i]-dipposx[j])+(dipposy[i]-dipposy[j])*(dipposy[i]-dipposy[j])+(dipposz[i]-dipposz[j])*(dipposz[i]-dipposz[j]));
            
            dipvecx=dipposx[i]-dipposx[j];
            dipvecy=dipposy[i]-dipposy[j];
            dipvecz=dipposz[i]-dipposz[j];
            
#ifdef DEBUG
            printf("Dipoles: %2d  %2d Difference vector: %5.3f %5.3f %5.3f  \nDifference vector length calc. 1: %5.3f\n",i,j,dipvecx,dipvecy,dipvecz,dipdist);
            printf("Dipole length: %2d = %5.3f, %2d = %5.3f \nDifference vector length calc. 2: %5.3f\n",i,length(dipposx[i],dipposy[i],dipposz[i]),j,length(dipposx[j],dipposy[j],dipposz[j]),length(dipposx[i]-dipposx[j],dipposy[i]-dipposy[j],dipposz[i]-dipposz[j]));
#endif
    
            if(dipdist<.2)printf ("\nDistances too small in interaction at %d and %d: %10.5lf\n",i,j,dipdist);
                
            // Write non-nearest neighbor elements of Hamiltonian with Transition Dipole Coupling (TDC) model:
            H1[i][j]=0;
            H1[i][j]+=( ( (vxd[i]*vxd[j])+(vyd[i]*vyd[j])+(vzd[i]*vzd[j]) ) / pow(dipdist,3) )-3*( ( ((dipvecx*vxd[i])+(dipvecy*vyd[i])+(dipvecz*vzd[i]))*((dipvecx*vxd[j])+(dipvecy*vyd[j])+(dipvecz*vzd[j]))) / (pow(dipdist,5)) );
    
    H1[i][j]=5033*H1[i][j];
    H1[j][i]=H1[i][j];
        }

}

free_f3tensor(xd,1,2,1,n,1,4);free_f3tensor(yd,1,2,1,n,1,4);free_f3tensor(zd,1,2,1,n,1,4);

}

'''





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
#x=matrix(1,n,1,5)
#vxd=vector(1,n)
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
            print i
        except (ValueError,IndexError):
            pass
    for i in atom_data:
        print i[6:9]
    n=len(atom_data)
    #In the x,y,z matrices defined below, we have the following:
    #   The first position is reserved for carbon atoms ('C')
    #   The second position is reserved for oxygen atoms ('O')
    #   The third position is reserved for nitrogen atoms ('N')
    #   The fourth position is reserved for the following: ('H','CN','HN')
    #   The fifth position is reserved for 
    x = [[0,0,0,0,0] for _ in range(n)]
    y = [[0,0,0,0,0] for _ in range(n)]
    z = [[0,0,0,0,0] for _ in range(n)]
    for i in range(0,len(atom_data)):
        if atom_data[i][2] == 'C':
            x[i][0] = float(atom_data[i][6])
            y[i][0] = float(atom_data[i][7])
            z[i][0] = float(atom_data[i][8])
        if atom_data[i][2] == 'O':
            x[i][1] = float(atom_data[i][6])
            y[i][1] = float(atom_data[i][7])
            z[i][1] = float(atom_data[i][8])
        if atom_data[i][2] == 'N':
            x[i][2] = float(atom_data[i][6])
            y[i][2] = float(atom_data[i][7])
            z[i][2] = float(atom_data[i][8])
        if atom_data[i][2] == 'N':
            x[i][3] = float(atom_data[i][6])
            y[i][3] = float(atom_data[i][7])
            z[i][3] = float(atom_data[i][8])
        if atom_data[i][2] == 'H' or atom_data[i][2] == 'CN' \
                    or atom_data[i][2] == 'HN':
            x[i][4] = float(atom_data[i][6])
            y[i][4] = float(atom_data[i][7])
            z[i][4] = float(atom_data[i][8])
    return atom_data,x,y,z,n


def ran1(idum):
    NTAB = 32
    IQ = 127773
    IA = 16807
    IR = 2836
    IM = 2147483647
    NDIV = 1+(IM-1)/NTAB
    AM = 1.0/IM
    RNMX = 1.0 - (1.2*10**-7)
    iy = 0
    iv = [0]*NTAB
    if idum<=0 or not iy:
        if -idum<1:
            idum = 1
        else:
            idum = -idum
        for j in range(NTAB+7,-1,-1):
            k = idum/IQ
            idum = IA*(idum-k*IQ)-IR*k
            if idum<0:
                idum += IM
            if j<NTAB:
                iv[j] = idum
        iy = iv[0]
    k = idum/IQ
    idum = IA*(idum-k*IQ)-IR*k
    if idum<0:
        idum+=IM
    j = int(iy/NDIV)
    iy = iv[j]
    iv[j] = idum
    if AM*iy>RNMX:
        return RNMX
    else:
        return AM*iy
    

def gasdev(idum):
    iset = 0
    gset = 0
    if iset == 0:
        t = 0
        while t>=0:
            v1 = 2.0*ran1(idum)-1.0
            v2 = 2.0*ran1(idum)-1.0
            rsq = v1*v1+v2*v2
            t = rsq
        fac = np.sqrt(-2.0*np.log(rsq)/rsq)
        gset = v1*fac
        iset = 1
        return v2*fac
    else:
        iset = 0
        return gset

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



def transition_charges(hamiltonian,vxd,vyd,vzd,x,y,z,n,istruct,protein):
    #CALCULATE THE DISPLACED ATOMS
    amplitude=(10**10)*np.sqrt(h/(8*math.pi*math.pi*nma_freq*c*nma_mass*amu))
    xd = [[[0]*n]*4]*2
    yd = [[[0]*n]*4]*2
    zd = [[[0]*n]*4]*2
         
    for i in range(0,n):
        vxCOrev = x_mode[0]-x_mode[1]
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
        
        for j in range(0,5):
            a_CO = (x_der[j]*vxCOrev+y_der[j]*vzCOrev)/length_3d(vxCOrev,vyCOrev,vzCOrev)/length_3d(vxCO,vyCO,vzCO)
            a_CN = (x_der[j]*vxCNrev+y_der[j]*vzCNrev)/length_3d(vxCNrev,vyCNrev,vzCNrev)/length_3d(vxCN,vyCN,vzCN)
            if j != 0 and a_CO>0:                
                xd[0][i][j] = x[i][j]+amplitude/2*(a_CO*vxCO+a_CN*vxCN)
                yd[0][i][j] = y[i][j]+amplitude/2*(a_CO*vyCO+a_CN*vyCN)
                zd[0][i][j] = z[i][j]+amplitude/2*(a_CO*vzCO+a_CN*vzCN)
                
                xd[1][i][j] = x[i][j]-amplitude/2*(a_CO*vxCO+a_CN*vxCN)
                yd[1][i][j] = y[i][j]-amplitude/2*(a_CO*vyCO+a_CN*vyCN)
                zd[1][i][j] = z[i][j]-amplitude/2*(a_CO*vzCO+a_CN*vzCN)
            else:
                xd[0][i][j] = x[i][j]-amplitude/2*(a_CO*vxCO+a_CN*vxCN)
                yd[0][i][j] = y[i][j]-amplitude/2*(a_CO*vyCO+a_CN*vyCN)
                zd[0][i][j] = z[i][j]-amplitude/2*(a_CO*vzCO+a_CN*vzCN)
                
                xd[1][i][j] = x[i][j]+amplitude/2*(a_CO*vxCO+a_CN*vxCN)
                yd[1][i][j] = y[i][j]+amplitude/2*(a_CO*vyCO+a_CN*vyCN)
                zd[1][i][j] = z[i][j]+amplitude/2*(a_CO*vzCO+a_CN*vzCN)
        #CALCULATE TRANSITION DIPOLES
        
    
        

'''
    for (i=1;i<=n;i++)
    {
        vxCOrev=xmode[1]-xmode[2];vyCOrev=ymode[1]-ymode[2];vzCOrev=zmode[1]-zmode[2]; ###
        vxCNrev=xmode[1]-xmode[3];vyCNrev=ymode[1]-ymode[3];vzCNrev=zmode[1]-zmode[3]; ###
        vxCO=x[i][1]-x[i][2];vyCO=y[i][1]-y[i][2];vzCO=z[i][1]-z[i][2]; ###
        vxCN=x[i][1]-x[i][3];vyCN=y[i][1]-y[i][3];vzCN=z[i][1]-z[i][3]; ###
        
        for (j=1;j<=4;j++)
        {
            aCO=(x_der[j]*vxCOrev+y_der[j]*vyCOrev+z_der[j]*vzCOrev)/length(vxCOrev,vyCOrev,vzCOrev)/length(vxCO,vyCO,vzCO);
            aCN=(x_der[j]*vxCNrev+y_der[j]*vyCNrev+z_der[j]*vzCNrev)/length(vxCNrev,vyCNrev,vzCNrev)/length(vxCN,vyCN,vzCN);
            if (j==1)
            {
                sign[i]=1;                                /*Test sign of Eigenmode by movement of C-atom)*/
                if (aCO>0) sign[i]=-1;
            }
            
            xd[1][i][j]=x[i][j]-sign[i]*Amplitude/2.*(aCO*vxCO+aCN*vxCN);
            yd[1][i][j]=y[i][j]-sign[i]*Amplitude/2.*(aCO*vyCO+aCN*vyCN);
            zd[1][i][j]=z[i][j]-sign[i]*Amplitude/2.*(aCO*vzCO+aCN*vzCN);
            xd[2][i][j]=x[i][j]+sign[i]*Amplitude/2.*(aCO*vxCO+aCN*vxCN);
            yd[2][i][j]=y[i][j]+sign[i]*Amplitude/2.*(aCO*vyCO+aCN*vyCN);
            zd[2][i][j]=z[i][j]+sign[i]*Amplitude/2.*(aCO*vzCO+aCN*vzCN);
        }
    }
'''

def parameterizedB3LYP(hamiltonian,x,y,z,n):
    ph = []
    ps = []
    for i in range(0,n):
        phi=dehydral(x[i][1],x[i][3],x[i][5],x[i+1][1],y[i][1],y[i][3],y[i][5],y[i+1][1],z[i][1],z[i][3],z[i][5],z[i+1][1])
        psi=dehydral(x[i][3],x[i][5],x[i+1][1],x[i+1][3],y[i][3],y[i][5],y[i+1][1],y[i+1][3],z[i][3],z[i][5],z[i+1][1],z[i+1][3])
        ph.append(phi)
        ps.append(psi)
        coupling = 0
        for k in range(0,7):
            for j in range(0,7):
                coupling += coupling_par[j+1+(k*7)]*math.cos(j*psi/180*math.pi)*math.cos(k*phi/180*math.pi)
        for k in range(0,6):
            for j in range(0,6):
                coupling += coupling_par[j+49+((k-1)*5)]*math.sin(j*psi/180*math.pi)*math.sin(k*phi/180*math.pi)
        hamiltonian[i][i+1] = coupling
        hamiltonian[i+1][i] = coupling #Building diagonal-adjacents to the hamiltonian matrix!
    






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
                
            

def HydrogenBondChimera(hamiltonian,n,subuoffset,subucumsum,prolist,istruct,protein):
    return

def HydrogenBondZanniSkinner(hamiltonian,n,subuoffset,subucumsum,prolist):
    return



def interaction(hamiltonian, vxd, vyd, vzd, x, y, z, n, name, \
                subuoffset, subucumsum, prolist, istruct, protein, \
                protocol = 'HydrogenBond'):
    '''Protocol may be HydrogenBond, HydrogenBondChimera, or HydrogenBondZanniSkinner.'''
    transition_charges(hamiltonian,vxd,vyd,vzd,x,y,z,n,istruct,protein)
    parameterizedB3LYP(hamiltonian,x,y,z,n)
    if protocol == 'HydrogenBond':
        HydrogenBond(hamiltonian,x,y,z,n,prolist)
    elif protocol == 'HydrogenBondChimera':
        #HydrogenBondChimera(hamiltonian,n,subuoffset,subucumsum,prolist,istruct,protein)
        raise NotImplementedError
    else:
        #HydrogenBondZanniSkinner(hamiltonian,n,subuoffset,subucumsum,prolist)
        raise NotImplementedError


def main_calculate(file_name):
    name = file_name.replace('.pdb','')
    atom_data,x,y,z,n = coordinate_read(file_name)
    if n == 0:
        raise ValueError('The PDB file appears to be empty. :(')
    hamiltonian = [[0]*n for _ in range(n)] #Instantiating the nxn matrix.
    #This will later be used for the Hamiltonian!
    for i in range(0,n):
        hamiltonian[i][i] += woff+w_inhomogeneous*gasdev(idum) #look into gasdev later...
    interaction(hamiltonian,vxd,vyd,vzd,x,y,z,n,name,subuoffset,subucumsum,prolist,istruct,protein)
        
main_calculate('/Users/mschwart/vsfg-bellerephon/2chb.pdb')
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
