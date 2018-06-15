import bellerephon_main as bm
import matplotlib.pyplot as plt
import os

#file_path = '/Users/mschwart/vsfg-bellerephon/test-16jan/'
#file_path = '/Users/mschwart/vsfg-bellerephon/alphas/confout/2bsk.1.15-deg/'

d = '/Users/mschwart/vsfg-bellerephon/alphas/confout/'
e = '/Users/mschwart/vsfg-bellerephon/betaproteins/confout/'

#subfolders = [f.path for f in os.scandir(d) if f.is_dir() ]
#print(subfolders)
t = os.walk(d)
s = os.walk(e)
#print([x[0] for x in s])
#full_list = [x[0] for x in s]


total_list = [['Wavelength (nm)']+[x * 0.5 for x in range(3200, 3601)]]


#print([x[0] for x in t])
pdb_list = []
trans_ham_list = []
eigenval_list = []
eigenvec_list = []

for file_path in [x[0] for x in s][1:]:
    protein_name = file_path.split('/')[-1]
    print('Now computing for...'+file_path)
    pdb_file = file_path+'/output.pdb'
    trans_ham_file = file_path+'/output_prepped_transdipmom_Hamm.txt'
    eigenval_file = file_path+'/output_prepped_eval.txt'
    eigenvec_file = file_path+'/output_prepped_evec.txt'
    iterval = 0
    if os.path.isfile(pdb_file):
        pdb_list.append(pdb_file)
    else:
        print 'PDB MISSING!'
        iterval+=1
        pdb_list.append(pdb_list[-1])

    if os.path.isfile(trans_ham_file):
        trans_ham_list.append(trans_ham_file)
    else:
        print 'TRANSITION HAMILTONIANS MISSING!'
        iterval+=1
        trans_ham_list.append(trans_ham_list[-1])

    if os.path.isfile(eigenval_file):
        eigenval_list.append(eigenval_file)
    else:
        print 'EIGENVALUES MISSING!'
        iterval+=1
        eigenval_list.append(eigenval_list[-1])
        
    if os.path.isfile(eigenvec_file):
        eigenvec_list.append(eigenvec_file)
    else:
        print 'EIGENVECTORS MISSING!'
        iterval+=1
        eigenvec_list.append(eigenvec_list[-1])
    
    if iterval < 2:
        OmegaVAL, VSFG_SPS, VSFG_SSP = bm.bellerephon_process(pdb_list[-1],trans_ham_list[-1],eigenval_list[-1],eigenvec_list[-1])
        plt.figure()
        plt.plot(OmegaVAL,VSFG_SPS)
        plt.plot(OmegaVAL,VSFG_SSP)
        plt.title('VSFG')
        plt.xlabel('Wavelength (nm)')
        plt.ylabel('Intensity (unitless)')
        plt.show()
        total_list.append([protein_name+'_SPS']+VSFG_SPS)
        total_list.append([protein_name+'_SSP']+VSFG_SSP)
        print 'Just finished '+protein_name
    else:
        print('FATAL ERROR! '+file_path)

import csv

total_list = zip(*total_list)


#with open('beta_proteins_vsfg.csv', "wb") as csv_file:
#    writer = csv.writer(csv_file, delimiter=',')
#    for line in total_list:
#        writer.writerow(line)
