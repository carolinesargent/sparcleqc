"""
Unit and regression test for the sparcle_qc package.
"""

# Import package, test suite, and other packages as needed
import sys
import shutil
import pytest
import os
import sparcle_qc

'''
@pytest.fixture(scope='session', autouse=True)
def cleanup_directory():
    #directory_path = 'test'

    # Remove the directory if it exists
    if os.path.exists(directory_path):
        shutil.rmtree(directory_path)

    yield  # This allows the tests to run
    # Cleanup after tests
    if os.path.exists(directory_path):
        print('removing files')
        shutil.rmtree(directory_path)
'''
       
def test_sparcle_qc_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "sparcle_qc" in sys.modules

def test_run_sapt_psi4_amber():
    inputs = {
    'input_filename': 'test1.in',
    'pdb_file': '4yff_cl_cleaned.pdb',
    'pre-capped': 'True',
    'cutoff': 8.5,
    'seed': 4247,
    'seed_file': '4yff_cl_cleaned.pdb',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'fisapt0',
    'fisapt_partition': 'True',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': '/theoryfs2/ds/ipberry/miniconda3/envs/emb_sapt/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run(user_options = inputs)
    true_dictionary = (392, 3903, 0.0, -0.0)
    assert output_dictionary == true_dictionary
def test_run_hf_psi4_amber():
    inputs = {
    'input_filename': 'test2.in',
    'pdb_file': '4yff_cl_cleaned.pdb',
    'pre-capped': 'True',
    'cutoff': 8.5,
    'seed': 4247,
    'seed_file': '4yff_cl_cleaned.pdb',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'hf',
    'fisapt_partition': 'True',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': '/theoryfs2/ds/ipberry/miniconda3/envs/emb_sapt/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run(user_options = inputs)
    true_dictionary = {'Complex': [392, 3903, 0.0, -0.0], 'Ligand': [19, 0, 0.0, 0.0], 'Protein': [373, 3903, 0.0, -0.0]}
    assert output_dictionary == true_dictionary
def test_run_sapt_psi4_charmm():
    inputs = {
    'input_filename': 'test3.in',
    'pdb_file': '3qxp.pdb',
    'cutoff': 5,
    'seed': 'ligand',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'fisapt0',
    'fisapt_partition': 'True',
    'basis_set': 'aug-cc-pv(D+d)z',
    'env_path': '/theoryfs2/ds/ipberry/miniconda3/envs/emb_sapt/',
    'charmm_rtf': 'top_all36_prot.rtf',
    'charmm_prm': 'par_all36m_prot.prm',
    'water_model': 'tip3p' ,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run(user_options =inputs)
    true_dictionary = (542, 4957, 1.0, 2.0)
    assert output_dictionary == true_dictionary
def test_run_hf_psi4_charmm():
    inputs = {
    'input_filename': 'test4.in',
    'pdb_file': '3qxp.pdb',
    'cutoff': 5,
    'seed': 'ligand',
    'charge_scheme': 'DZ3',
    'ligand_charge': 0,
    'method': 'hf',
    'fisapt_partition': 'True',
    'basis_set': 'aug-cc-pv(D+d)z',
    'env_path': '/theoryfs2/ds/ipberry/miniconda3/envs/emb_sapt/',
    'charmm_rtf': 'top_all36_prot.rtf',
    'charmm_prm': 'par_all36m_prot.prm',
    'water_model': 'tip3p' ,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run(user_options = inputs)
    true_dictionary = {'Complex': [542, 4845, 1.0, 2.0], 'Ligand': [41, 0, 0.0, 0.0], 'Protein': [501, 4845, 1.0, 2.0]}
    assert output_dictionary == true_dictionary

def test_run_sapt_qchem_amber():
    inputs = {
    'input_filename': 'test5.in',
    'pdb_file': '4yff_cl_cleaned.pdb',
    'pre-capped': 'True',
    'cutoff': 8.5,
    'seed': 4247,
    'seed_file': '4yff_cl_cleaned.pdb',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'sapt0',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': '/theoryfs2/ds/ipberry/miniconda3/envs/emb_sapt/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'q-chem',
    'mem': '60 GB',
    'nthreads': 2}
    output_dictionary= sparcle_qc.run(user_options = inputs)
    true_dictionary = (392, 3903, 0.0, -0.0)
    assert output_dictionary == true_dictionary
def test_run_hf_qchem_amber():
    inputs = {
    'input_filename': 'test6.in',
    'pdb_file': '4yff_cl_cleaned.pdb',
    'pre-capped': 'True',
    'cutoff': 8.5,
    'seed': 4247,
    'seed_file': '4yff_cl_cleaned.pdb',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': '/theoryfs2/ds/ipberry/miniconda3/envs/emb_sapt/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'q-chem',
    'mem': '60 GB',
    'nthreads': 2,
    'cp': 'false'}
    output_dictionary= sparcle_qc.run(user_options = inputs)
    true_dictionary = {'Complex': [392, 3903, 0.0, -0.0], 'Ligand': [19, 0, 0.0, 0.0], 'Protein': [373, 3903, 0.0, -0.0]}
    assert output_dictionary == true_dictionary
def test_run_sapt_qchem_charmm():
    inputs = {
    'input_filename': 'test7.in',
    'pdb_file': '3qxp.pdb',
    'cutoff': 5,
    'seed': 'ligand',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'sapt0',
    'basis_set': 'aug-cc-pv(D+d)z',
    'env_path': '/theoryfs2/ds/ipberry/miniconda3/envs/emb_sapt/',
    'charmm_rtf': 'top_all36_prot.rtf',
    'charmm_prm': 'par_all36m_prot.prm',
    'water_model': 'tip3p' ,
    'software': 'q-chem',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run(user_options =inputs)
    true_dictionary = (542, 4957, 1.0, 2.0)
    assert output_dictionary == true_dictionary
def test_run_hf_qchem_charmm():
    inputs = {
    'input_filename': 'test8.in',
    'pdb_file': '3qxp.pdb',
    'cutoff': 5,
    'seed': 'ligand',
    'charge_scheme': 'DZ3',
    'ligand_charge': 0,
    'method': 'hf',
    'fisapt_partition': 'True',
    'basis_set': 'aug-cc-pv(D+d)z',
    'env_path': '/theoryfs2/ds/ipberry/miniconda3/envs/emb_sapt/',
    'charmm_rtf': 'top_all36_prot.rtf',
    'charmm_prm': 'par_all36m_prot.prm',
    'water_model': 'tip3p' ,
    'software': 'q-chem',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run(user_options = inputs)
    true_dictionary = {'Complex': [542, 4845, 1.0, 2.0], 'Ligand': [41, 0, 0.0, 0.0], 'Protein': [501, 4845, 1.0, 2.0]}
    assert output_dictionary == true_dictionary
def test_run_hf_nwchem_amber():
    inputs = {
    'input_filename': 'test9.in',
    'pdb_file': '4yff_cl_cleaned.pdb',
    'pre-capped': 'True',
    'cutoff': 8.5,
    'seed': 4247,
    'seed_file': '4yff_cl_cleaned.pdb',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': '/theoryfs2/ds/ipberry/miniconda3/envs/emb_sapt/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'nwchem',
    'nwchem_scratch': '/scratch/ipberry',
    'nwchem_perm': '/scratch/ipberry',
    'mem': '60 GB',
    'nthreads': 2,
    'cp': 'false'}
    output_dictionary= sparcle_qc.run(user_options = inputs)
    true_dictionary = {'Complex': [392, 3903, 0.0, -0.0], 'Ligand': [19, 0, 0.0, 0.0], 'Protein': [373, 3903, 0.0, -0.0]}
    assert output_dictionary == true_dictionary
def test_run_hf_nwchem_charmm():
    inputs = {
    'input_filename': 'test10.in',
    'pdb_file': '3qxp.pdb',
    'cutoff': 5,
    'seed': 'ligand',
    'charge_scheme': 'DZ3',
    'ligand_charge': 0,
    'method': 'hf',
    'fisapt_partition': 'True',
    'basis_set': 'aug-cc-pv(D+d)z',
    'env_path': '/theoryfs2/ds/ipberry/miniconda3/envs/emb_sapt/',
    'charmm_rtf': 'top_all36_prot.rtf',
    'charmm_prm': 'par_all36m_prot.prm',
    'water_model': 'tip3p' ,
    'software': 'nwchem',
    'nwchem_scratch': '/scratch/ipberry',
    'nwchem_perm': '/scratch/ipberry',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run(user_options = inputs)
    true_dictionary = {'Complex': [542, 4845, 1.0, 2.0], 'Ligand': [41, 0, 0.0, 0.0], 'Protein': [501, 4845, 1.0, 2.0]}
    assert output_dictionary == true_dictionary
