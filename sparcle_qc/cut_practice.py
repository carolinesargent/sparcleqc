from cut_protein import run_cut_protein 
import os

os.chdir('inputs')
run_cut_protein('cx_autocap_fixed.pdb', 'ligand', 5)
