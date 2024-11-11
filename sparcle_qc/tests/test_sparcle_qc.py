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
    'env_path': f'{os.getcwd()}/',
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
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
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
    'charmm_rtf': 'top_all36_prot.rtf',
    'charmm_prm': 'par_all36m_prot.prm',
    'water_model': 'tip3p' ,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
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
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'q-chem',
    'mem': '60 GB',
    'nthreads': 2}
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
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
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'q-chem',
    'mem': '60 GB',
    'nthreads': 2,
    'cp': 'false'}
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
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
    'charmm_rtf': 'top_all36_prot.rtf',
    'charmm_prm': 'par_all36m_prot.prm',
    'water_model': 'tip3p' ,
    'software': 'q-chem',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run_sparcle(user_options =inputs)
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
    'charmm_rtf': 'top_all36_prot.rtf',
    'charmm_prm': 'par_all36m_prot.prm',
    'water_model': 'tip3p' ,
    'software': 'q-chem',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
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
    'env_path': f'{os.getcwd()}/',
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
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
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
    'charmm_rtf': 'top_all36_prot.rtf',
    'charmm_prm': 'par_all36m_prot.prm',
    'water_model': 'tip3p' ,
    'software': 'nwchem',
    'nwchem_scratch': '/scratch/ipberry',
    'nwchem_perm': '/scratch/ipberry',
    'mem': '100 GB',
    'nthreads': 16,
    'cp': 'true'}
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = {'Complex': [542, 4845, 1.0, 2.0], 'Ligand': [41, 0, 0.0, 0.0], 'Protein': [501, 4845, 1.0, 2.0]}
    assert output_dictionary == true_dictionary
def test_run_dz1():
    inputs = {
    'input_filename': 'test11.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 5,
    'seed': 'ligand',
    'charge_scheme': 'DZ1',
    'ligand_charge': 0,
    'method': 'sapt0',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }
    
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = (586, 4330, 4.0, 2.0)
    assert output_dictionary == true_dictionary
def test_run_dz2():
    inputs = {
    'input_filename': 'test12.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 5,
    'seed': 'ligand',
    'charge_scheme': 'DZ2',
    'ligand_charge': 0,
    'method': 'sapt0',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }
    
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = (586, 4290, 4.0, 2.0)
    assert output_dictionary == true_dictionary
def test_run_z1():
    inputs = {
    'input_filename': 'test13.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'Z1',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }
    
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = {'Complex': [721, 4207, 2.0, -0.3], 'Ligand': [41, 0, 0.0, 0.0], 'Protein': [680, 4207, 2.0, -0.3]}
    assert output_dictionary == true_dictionary
def test_run_z2():
    inputs = {
    'input_filename': 'test14.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'Z2',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }
    
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = {'Complex': [721, 4207, 2.0, 6.97], 'Ligand': [41, 0, 0.0, 0.0], 'Protein': [680, 4207, 2.0, 6.97]}
    assert output_dictionary == true_dictionary
def test_run_z3():
    inputs = {
    'input_filename': 'test15.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'Z3',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }
    
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = {'Complex': [721, 4207, 2.0, 4.97], 'Ligand': [41, 0, 0.0, 0.0], 'Protein': [680, 4207, 2.0, 4.97]} 
    assert output_dictionary == true_dictionary
def test_run_brcd():
    inputs = {
    'input_filename': 'test16.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'BRCD',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }
    
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = {'Complex': [721, 4225, 2.0, 3.0], 'Ligand': [41, 0, 0.0, 0.0], 'Protein': [680, 4225, 2.0, 3.0]}
    assert output_dictionary == true_dictionary
def test_run_brc2():
    inputs = {
    'input_filename': 'test17.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'BRC2',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }
    
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = {'Complex': [721, 4195, 2.0, 3.0], 'Ligand': [41, 0, 0.0, 0.0], 'Protein': [680, 4195, 2.0, 3.0]}
    assert output_dictionary == true_dictionary
def test_run_see():
    inputs = {
    'input_filename': 'test18.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'SEE',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }
    
    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = {'Complex': [721, 78, 2.0, -1.81], 'Ligand': [41, 0, 0.0, 0.0], 'Protein': [680, 78, 2.0, -1.81]}
    assert output_dictionary == true_dictionary
def test_exit():
    exits = []
    true_exits = []
    inputs = [{
    'input_filename': 'test20.in',
    'pdb_file': 'doesnt_exit.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    },{'input_filename': 'test21.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'ABC',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    },{'input_filename': 'test22.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 'six',
    'seed': 'ligand',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    },{'input_filename': 'test23.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    },{'input_filename': 'test24.in',
    'pdb_file': '3qxp.pdb',
    'cutoff': '6',
    'seed': 'ligand',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'charmm_rtf': 'top_all36_prot.rtf',
    'env_path': f'{os.getcwd()}/'
    'water_model': 'tip3p' ,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }]
    for input_test in inputs:
        with pytest.raises(SystemExit) as e: 
            output_dictionary= sparcle_qc.run_sparcle(user_options = input_test)
        exits.append(e.type)
        true_exits.append(SystemExit)
        
    assert exits == true_exits 

def test_run_convert():
    inputs = {
    'input_filename': 'test19.in',
    'pdb_file': '3QU0_templated_from_3QXP_amber.pdb',
    'template_path': 'reference_convert/cx_autocap_fixed.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'BRC',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }

    output_dictionary= sparcle_qc.run_sparcle(user_options = inputs)
    true_dictionary = {'Complex': [718, 4225, 2.0, 3.0], 'Ligand': [38, 0, 0.0, 0.0], 'Protein': [680, 4225, 2.0, 3.0]}
    assert output_dictionary == true_dictionary

