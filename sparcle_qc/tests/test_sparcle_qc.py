"""
Unit and regression test for the sparcle_qc package.
"""

# Import package, test suite, and other packages as needed
import sys
import shutil
import pytest
import os
import sparcle_qc

def delete_test_files():
    for file in os.listdir('.'):
        if file.startswith('test') and file.endswith('.in'):
            os.remove(file)
            print(f"Removed file: {file}")

def delete_test_dirs():
    for dir_name in os.listdir('.'):
        if dir_name.startswith('test') and os.path.isdir(dir_name):
            shutil.rmtree(dir_name)
            print(f"Removed directory: {dir_name}")

@pytest.fixture(scope='session', autouse=True)
def cleanup_directory():
    yield  # This allows the tests to run
    # Cleanup after tests
    delete_test_files()
    delete_test_dirs()


def test_sparcle_qc_imported():
    """Sample test, will always pass so long as import statement worked."""
    assert "sparcle_qc" in sys.modules

def test_run_sapt_psi4_amber():
    print(os.getcwd())
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
    'env_path': '/home/runner/work/sparcle_qc/sparcle_qc/sparcle_qc/tests/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = (392, 3903, 0.0, -0.0)
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
    'charmm_rtf': 'top_all36_prot.rtf',
    'charmm_prm': 'par_all36m_prot.prm',
    'water_model': 'tip3p' ,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run_sparcle(user_options =inputs)
    true_dictionary = (542, 4957, 1.0, 2.0)
    assert output_dictionary == true_dictionary

