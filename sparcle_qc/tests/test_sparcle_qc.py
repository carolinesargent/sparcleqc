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
@pytest.fixture(autouse = True)
def dictionary_inputs():
    inputs = {
    'input_filename': 'dictionary_test.in',
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
    yield inputs
       
def test_sparcle_qc_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "sparcle_qc" in sys.modules

def test_run_sapt_psi4_amber(dictionary_inputs):
    dictionary_inputs['input_filename'] = 'test1.in'
    output_dictionary= sparcle_qc.run(user_options = dictionary_inputs)
    true_dictionary = (392, 3903, 0.0, -0.0)
    assert output_dictionary == true_dictionary
def test_run_hf_psi4_amber(dictionary_inputs):
    dictionary_inputs['method'] = 'hf'
    dictionary_inputs['input_filename'] = 'test2.in'
    output_dictionary= sparcle_qc.run(user_options = dictionary_inputs)
    true_dictionary = {'Complex': [392, 3903, 0.0, -0.0], 'Ligand': [19, 0, 0.0, 0.0], 'Protein': [373, 3903, 0.0, -0.0]}
    assert output_dictionary == true_dictionary
def test_run_sapt_psi4_charmm(dictionary_inputs):
    dictionary_inputs['input_filename'] = 'test3.in'
    dictionary_inputs['pdb_file'] = '3qxp.pdb'
    dictionary_inputs['seed'] = 'ligand'
    dictionary_inputs['cutoff'] = 5
    dictionary_inputs.pop('pre-capped')
    dictionary_inputs.pop('seed_file')
    dictionary_inputs.pop('amber_ff')
    dictionary_inputs['charmm_rtf'] = 'top_all36_prot.rtf'
    dictionary_inputs['charmm_prm'] = 'par_all36m_prot.prm'
    dictionary_inputs['water_model'] = 'tip3p'
    dictionary_inputs.pop('o_charge')
    dictionary_inputs.pop('h_charge')
    dictionary_inputs.pop('ep_charge')
    output_dictionary= sparcle_qc.run(user_options = dictionary_inputs)
    true_dictionary = (542, 4957, 1.0, 2.0)
    assert output_dictionary == true_dictionary

