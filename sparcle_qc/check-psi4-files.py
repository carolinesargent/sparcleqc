import numpy as np
import os

psi4files =[x for x in os.listdir() if '-psi4.py' in x]
#psi4files = ['cut-Z1-psi4.py']

for psi4file in psi4files:
    with open(psi4file, 'r') as pfile:
        lines = pfile.readlines()
        extern_idx = []
        dimer_idx = []
        prot_idx = []
        for n,l in enumerate(lines):
            if 'Chargefield' in l:
                extern_idx.append(n)
            elif 'dimer' in l:
                dimer_idx.append(n)
            elif '--' in l:
                prot_idx.append(n)
            elif 'unit' in l:
                dimer_idx.append(n)
        array = lines[extern_idx[0]+1:extern_idx[1]]
        charge = float(array[0].split(',')[0])
        num_atoms = dimer_idx[1] - dimer_idx[0] - 4 + 1
        num_prot_atoms = dimer_idx[1] - prot_idx[0] - 2
        for l in array[1:]:
            charge += float(l.split(',')[1])
            num_atoms += 1
        print(f'{psi4file}: {charge:.2f}, {num_atoms}')
        print(f'Number of protein atoms: {num_prot_atoms}')


