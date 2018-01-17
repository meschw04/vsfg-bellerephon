#HAMILTONIAN CALCULATIONS!!!
import math
import numpy as np

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

'''
nma.par

nma_freq = 1746.0686
nma_mass = 9.8478
 0.483   0.160  -0.001   0.603948
 0.385   1.385  -0.000  -0.534156
-0.616  -0.650  -0.010  -0.515416
-0.492  -1.651  -0.020   0.390304
-0.03    0.73    0.00    0.007716
 0.04   -0.43    0.00    0.018198
-0.03   -0.07    0.00   -0.026049
-0.10   -0.10    0.00    0.000135
'''



'''
coupling.par

2.98049155409105
-6.599810192233803
-0.30853721377655763
0.08082590050798008
0.04740097894564941
0.008048225450833241
0.0015733734467524448
-0.9658675048030658
-5.22997717307316
-0.4018105791392881
-0.017339459064999913
0.008386549055078336
-0.050489074051387244
0.006789470425119076
-0.3126564488007089
-3.7797746273994806
-0.11414803857970511
-0.00017611006675912795
0.12579542585907855
0.011124863033873535
0.013850235703546394
-0.029503792846472005
-0.5059330170060446
0.19249211456707013
0.04314979965266982
0.15582397653156857
0.00007122142283001677
0.03310964759175535
-0.05620365560427052
-0.09323618884490228
0.07271537246962877
-0.006111394586803572
0.1140144332728223
-0.030650858533796854
0.010434624767047847
-0.006201344264232881
-0.07180433223027921
0.040607634844420835
-0.0979541787497221
0.11934199604608554
-0.012207576981502277
-0.018422318034652232
0.01883823305948914
0.0424046662659559
-0.03887914582205208
-0.1022335472962132
0.07300790795664054
-0.08500015077795386
-0.04615152341898034
6.610403410038493
0.590712804631773
-0.24362981965778352
-0.08002649779702173
0.019711383777822555
4.993250818063718
0.17452844187043454
-0.16960105630340355
-0.06764409458606896
-0.013064547947709688
0.27881995936872217
-0.3207748042878569
-0.03773019256433872
-0.10820787738659833
-0.05028414650455027
0.02492705580043824
0.01010521093108222
0.021042805555903196
-0.018502096344155176
-0.05345701390359108
0.06185935268126845
-0.01716502455463741
0.050050157280630725
-0.0820698925785323
-0.04129646850913813
'''



'''
/* Need this one! */
void transition_charges(float **H1,float *vxd,float *vyd,float *vzd,float **x,float **y,float **z,int n, int istruct, const char *protein)
{
    char temp[100];
    float xmode[5],ymode[5],zmode[5],qmode[5],Freq,Mass;
    float x_der[5],y_der[5],z_der[5],q_der[5];
    float vxCOrev,vyCOrev,vzCOrev,vxCO,vyCO,vzCO;
    float *dipposx,*dipposy,*dipposz;
    float vxCNrev,vyCNrev,vzCNrev,vxCN,vyCN,vzCN;
    float aCO,aCN,Amplitude,xtmp,q,sig,q1,q2,q3,q4,r1,r2,r3,r4;
    float ***xd, ***yd, ***zd;
    FILE *file,*transdipmomfile;
    int i,i1,i2,j,j1,j2,l,sign[nmax];
    
    xd=f3tensor(1,2,1,n,1,4);yd=f3tensor(1,2,1,n,1,4);zd=f3tensor(1,2,1,n,1,4);
    dipposx=vector(1,n);dipposy=vector(1,n);dipposz=vector(1,n);

    
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





'''
### SEE LINE 123!!!
float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
    long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
    float ***t;
    
    /* allocate pointers to pointers to rows */
    t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
    if (!t) nrerror("allocation failure 1 in f3tensor()");
    t += NR_END;
    t -= nrl;
    
    /* allocate pointers to rows and set pointers to them */
    t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
    if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
    t[nrl] += NR_END;
    t[nrl] -= ncl;
    
    /* allocate rows and set pointers to them */
    t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
    if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
    t[nrl][ncl] += NR_END;
    t[nrl][ncl] -= ndl;
    
    for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
    for(i=nrl+1;i<=nrh;i++) {
        t[i]=t[i-1]+ncol;
        t[i][ncl]=t[i-1][ncl]+ncol*ndep;
        for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
    }
    
    /* return pointer to array of pointers to rows */
    return t;
}

'''



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


def transition_charges(hamiltonian,vxd,vyd,vzd,x,y,z,n,istruct,protein):
#    f3tensor(1,2,1,n,1,4)
    amplitude=1e10*np.sqrt(h/(8*math.pi*math.pi*nma_freq*c*nma_mass*amu));

def parameterizedB3LYP(hamiltonian,x,y,z,n,subuchange):
    return

def HydrogenBond(hamiltonian,x,y,z,n,prolist):
    return    

def HydrogenBondChimera(hamiltonian,n,subuoffset,subucumsum,prolist,istruct,protein):
    return

def HydrogenBondZanniSkinner(hamiltonian,n,subuoffset,subucumsum,prolist):
    return


'''
void interaction(float **H1,float *vxd,float *vyd,float *vzd,float **x,float **y,float **z,int n,char *name, int *subuchange, int *subuoffset, int *subucumsum, int *prolist, int istruct, const char *protein)
{
    int i,j;
    
    for(i=1;i<=n;i++) for(j=1;j<=n;j++) H1[i][j]=0;
    
    transition_charges(H1,vxd,vyd,vzd,x,y,z,n,istruct,protein); /* Use transition charges from ab-initio calculations; this also calculates dipoles*/
    
    parametrizedB3LYP(H1,x,y,z,n,subuchange);	/* Use a pre-parametrized B3LYP calculation for nearest neighbors*/ 
    
    Hydrogenbond(H1,x,y,z,n,prolist);                   /* Search for hydrogen bonds and calculate frequency shifts directly from PDB-file (thus neglecting HBs to/from solvating water molecules and side chains) */
    
//  HydrogenbondChimera(H1,n,subuoffset,subucumsum,prolist,istruct,protein); /* Calculate frequency shifts from hydrogen bonds determined in Chimera */
    
//  HydrogenbondZanniSkinner(H1,n,subuoffset,subucumsum,prolist); /* Calculate frequency shifts from hydrogen bonds based on shifts in the electric field surrounding the peptide bond */
}
'''

def interaction(hamiltonian, vxd, vyd, vzd, x, y, z, n, name, subuchange, \
                subuoffset, subucumsum, prolist, istruct, protein, \
                protocol = 'HydrogenBond'):
    '''Protocol may be HydrogenBond, HydrogenBondChimera, or HydrogenBondZanniSkinner.'''
    transition_charges(hamiltonian,vxd,vyd,vzd,x,y,z,n,istruct,protein)
    parameterizedB3LYP(hamiltonian,x,y,z,n,subuchange)
    if protocol == 'HydrogenBond':
        HydrogenBond(hamiltonian,x,y,z,n,prolist)
    elif protocol == 'HydrogenBondChimera':
        HydrogenBondChimera(hamiltonian,n,subuoffset,subucumsum,prolist,istruct,protein)
    else:
        HydrogenBondZanniSkinner(hamiltonian,n,subuoffset,subucumsum,prolist)
