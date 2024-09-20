import sys
import shutil
import subprocess
import os
import threading
import itertools
import time
from typing import Dict

from .amber_prep import write_cpptraj, write_cpptraj_skip_autocap, write_tleap, autocap, skip_autocap, fix_numbers_amber
from .charmm_prep import dictionary_nocut, psf_to_mol2, combine_charmm, fix_numbers_charmm
from .pdb_prep import check_df_charges, check_resi_charges, convert_atom_id 
from .combine_data import combine_data
from .cut_protein import run_cut_protein 
from .convert_dict import convert_dictionary 
from .move_M3s import move_m3s
from .df_make_psi4 import write_QM, check_QM_file, write_input 
from .cap import run_cap
from .make_partition import partition 

stop_flashing = threading.Event()

def flashing_sparkle() -> None:
    """
    creates the flashing sparkle on the terminal

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    sparkle_emoji = "✨"
    cycle = itertools.cycle([sparkle_emoji, ' '])
    
    # Loop until the stop_flashing event is set
    while not stop_flashing.is_set():
        sys.stdout.write(next(cycle) + '\r')
        sys.stdout.flush()
        time.sleep(0.25)

def input_parser(filename:str) -> Dict:
    """
    parses the specified input file into a dictionary of keywords
    ensures that the specified inputs are valid options
    ensures that the necessary inputs are specified

    Parameters
    ----------
    filename: str
        path to input file

    Returns
    -------
    keywords: Dict
        dictionary of specified inputs
    """
    keywords = {}
    with open(filename, 'r') as f:
        for line in f:
            line = line.split('#')[0]

            if line not in ["", '\n']:
                split_line = line.split(':')
                key_word = split_line[0].strip()
                try:
                    value = split_line[1].strip()
                except IndexError:
                    print('Error: Invalid input file. Check that your options and values are separated by a colon\n')
                    sys.exit()
                if key_word == 'pdb_file':
                    if os.path.isfile(value) == False:
                        print('Error: Invalid input file. Path to PDB does not exist')
                        sys.exit()
                if key_word == 'pre-capped':
                    value = value.lower()
                    if value!= 'true' and value!='false':
                        print("Error: Invalid input file. Pre-capped is not true or false")
                        sys.exit()
                if key_word == 'cutoff':
                    try:
                        cutoff = float(value)
                    except ValueError:
                        print("Error: Invalid input file. Cutoff radius is not a numerical value")
                        sys.exit()
                if key_word =='seed':
                    try:
                        seed = int(value)
                    except ValueError:
                        if value != 'ligand':
                            print("Error: Invalid input file. Seed is not a numerical value or 'ligand'")
                            sys.exit()
                if key_word == 'seed_file':
                    if os.path.isfile(value) == False:
                        print('Error: Invalid input file. Path to seed_file PDB does not exist')
                        sys.exit()
                if key_word =='charge_scheme':
                    if value not in ['Z1', 'Z2','Z3','DZ1','DZ2','DZ3','BRC','BRCD','BRC2', 'SEE']:
                        print('Error: Invalid input file. Specified charge scheme is not a supported option')
                        sys.exit()
                if key_word =='ligand_charge':
                    try:
                        charge = int(value)
                    except ValueError:
                        print("Error: Invalid input file. Ligand charge is not an integer")
                        sys.exit()
                if key_word == 'fisapt_partition':
                    value = value.lower()
                    if value!= 'true' and value!='false':
                        print("Error: Invalid input file. fisapt_partition is not true or false")
                        sys.exit()
                if key_word == 'amber_ff':
                    #potentially check if the ff is in amber
                    pass
                if key_word =='env_path':
                    if os.path.exists(value) == False:
                        print('Error: Invalid input file. Path to environment does not exist')
                        sys.exit()
                if key_word == 'charmm_rtf':
                    if os.path.isfile(value) == False:
                        print('Error: Invalid input file. Path to CHARMM topology (.rtf) does not exist')
                        sys.exit()
                if key_word =='charmm_prm':
                    if os.path.isfile(value) == False:
                        print('Error: Invalid input file. Path to CHARMM parameters (.prm) does not exist')
                        sys.exit()
                if key_word =='water_model':
                    #potentially check later with tleap?
                    pass
                if key_word =='o_charge':
                    try:
                        charge = float(value)
                    except ValueError:
                        print("Error: Invalid input file. Custom water charge for oxygen is not a numerical value")
                        sys.exit()
                if key_word =='h_charge':
                    try:
                        charge = float(value)
                    except ValueError:
                        print("Error: Invalid input file. Custom water charge for hydrogen is not a numerical value")
                        sys.exit()
                if key_word =='ep_charge':
                    try:
                        charge = float(value)
                    except ValueError:
                        print("Error: Invalid input file. Custom extra water point charge is not a numerical value")
                        sys.exit()

                if key_word == 'template_path':
                    if os.path.isfile(value) == False:
                        print('Error: Invalid input file. Path to template PDB does not exist')
                        sys.exit()
                if key_word == 'do_fsapt':
                    value = value.lower()
                    if value != 'true' and value != 'false':
                        print("Error: Invalid input file. do_fsapt is not true or false")
                        sys.exit()
                keywords[key_word]=value
            
    if 'cutoff' not in keywords.keys():
        print('Error: Invalid input file. No cutoff is provided')
        sys.exit()
    if 'seed' not in keywords.keys():
        print('Error: Invalid input file. No seed is provided')
        sys.exit()
    if keywords['seed'] != 'ligand':
        if 'seed_file' not in keywords.keys():
            print('Error: Invalid input file. No seed file is provided for single atom seed')
            sys.exit()
    if 'charge_scheme' not in keywords.keys():
        print('Error: Invalid input file. No cutoff is provided')
        sys.exit()
    if 'method' not in keywords.keys():
        print('Error: Invalid input file. No cutoff is provided')
        sys.exit()
    if 'basis_set' not in keywords.keys():
        print('Error: Invalid input file. No cutoff is provided')
        sys.exit()
    if 'amber_ff' not in keywords.keys() and 'charmm_rtf' not in keywords.keys():
        print('Error: Invalid input file. CHARMM or Amber FF is not provided')
        sys.exit()
    if 'amber_ff' in keywords.keys():
        if 'env_path' not in keywords.keys():
            print('Error: Invalid input file. Path to enviornment is not provided')
            sys.exit()
        if 'water_model' not in keywords.keys():
            print('Error: Invalid input file. Water model is not provided')
            sys.exit()
    if 'amber_ff' in keywords.keys() and ('charmm_rtf' in keywords.keys() or 'charmm_prm' in keywords.keys()):
        print('Error: Invalid input file. CHARMM and Amber forcefields are provided')
        sys.exit()
    if 'charmm_rtf' in keywords.keys() or 'charmm_prm' in keywords.keys():
        try:
            charmm_prm = keywords['charmm_prm']
            charmm_rtf = keywords['charmm_rtf']
        except KeyError:
            print('Error: Invalid input file. CHARMM topology AND parameters are not provided')
            sys.exit()
    if 'h_charge' in keywords.keys() or 'o_charge' in keywords.keys():
        try:
            o_charge = keywords['o_charge']
            h_charge = keywords['h_charge']
        except KeyError:
            print('Error: Invalid input file. Both Oxygen and Hydrogen charges are not provided for water')
            sys.exit()
    if 'do_fsapt' not in keywords.keys():
        keywords['do_fsapt'] = 'false'
    if 'ep_charge' in keywords.keys():
        try:
            o_charge = keywords['o_charge']
            h_charge = keywords['h_charge']
        except KeyError:
            print('Error: Invalid input file. Oxygen and/or Hydrogen charges are not provided for water')
            sys.exit()

    print(f"\u2728Sparcle-QC is sparkling\u2728\nBeginning file preparation for an embedded QM calculation of {keywords['pdb_file']} ")
    
    return keywords


def run(input_file) -> None:
    """
    given an input file, parses the specified parameters into a dictionary of keywords and runs the necessary sparcle_qc functions to create an input file to be run with a quantum chemistry software by preparing pdbs, obtaining mol2s, carving out the QM region, capping the Q1 bonds with hydrogens, and finally writing an input file

    Parameters
    ----------
    input_file: str
        path to input file

    Returns
    -------
    None
    """
    
    #starting sparkle on command line
    flashing_thread = threading.Thread(target=flashing_sparkle)
    flashing_thread.start()
    try:
        #parsing input file into dictionary
        keywords = input_parser(input_file)
        #creating new directory for the created files, changing working directories, and copying necessary files into the new directory
        new_dir = input_file[:-3]
        os.mkdir(new_dir)
        if 'charmm_rtf' in keywords:
            shutil.copy(keywords['pdb_file'], new_dir)
            shutil.copy(keywords['pdb_file'].replace('pdb', 'psf'), new_dir)
            shutil.copy('ligand.pdb', new_dir)
        else:
            shutil.copy(keywords['pdb_file'], new_dir)
        shutil.move(input_file, new_dir)
        os.chdir(new_dir)
        output = open(f'{new_dir}.out', 'w')
        output.write('----------------------------------------------------------------------------------------------------\n')
        output.write('''                                                                                    QC 
                                                                                    /
                                                                                   /  *
                   ____                       _             ___   ____      /\u203E\u203E\u203E\u203E\u203E\u203E\\    *
                  / ___| _ __   __ _ _ __ ___| | ___       / _ \\ / ___|    /  \u00B7\u00B7\u00B7\u00B7  \\ *       
                  \\___ \\| '_ \\ / _` | '__/ __| |/ _ \\_____| | | | |       /  \u00B7    \u00B7  \\        
                   ___) | |_) | (_| | | | (__| |  __/_____| |_| | |___    \\  \u00B7    \u00B7  /  
                  |____/| .__/ \\__,_|_|  \\___|_|\\___|      \\__\\_\\\\____| *  \\  \u00B7\u00B7\u00B7\u00B7  /  
                        |_|                                               * \\______/ 
                                                                        *    \n''')
        output.write('----------------------------------------------------------------------------------------------------\n\n\n')
        #if forcefield is amber, writing and running cpptraj files and dealing with capping residues
        if 'amber_ff' in keywords: 
            if 'pre-capped' in keywords:
                if keywords['pre-capped'] == 'true':
                    write_cpptraj_skip_autocap(keywords['pdb_file'])
                else:
                    write_cpptraj(keywords['pdb_file'])
            else:
                write_cpptraj(keywords['pdb_file'])
            result = subprocess.run(['cpptraj -i cpptraj.in'], text = True, shell = True, capture_output = True)
            output.write('----------------------------------------------------------------------------------------------------\n')
            output.write('cpptraj'.center(100)+'\n')
            output.write('----------------------------------------------------------------------------------------------------\n')
            output.write(result.stdout)
            if 'pre-capped' in keywords:
                if keywords['pre-capped'] == 'true':
                    skip_autocap('cx_autocap.pdb')
                else:
                    autocap('uncapped.pdb')
            else:
                autocap('uncapped.pdb')
            fix_numbers_amber('prot_autocap.pdb')
            fix_numbers_amber('cx_autocap.pdb')
            os.remove('prot_autocap.pdb')
        #otherwise, the forcefield is charmm and combining protein and ligand into a complex pdb
        else:
            combine_charmm(keywords['pdb_file'])
            fix_numbers_charmm('cx_autocap.pdb')
        os.remove('cx_autocap.pdb')

        #obtaining seed containing information of which group to grow the QM region from
        if keywords['seed'] =='ligand':
            seed = 'ligand'
        else:
            seed = convert_atom_id(keywords['seed'], keywords['seed_file'])
        
        #if forcefield is amber, writing and running tleap
        if 'amber_ff' in keywords:
            write_tleap(keywords['amber_ff'], keywords['water_model'])
            result = subprocess.run(['tleap -f tleap.in'], text = True, shell = True, capture_output = True)
            output.write('----------------------------------------------------------------------------------------------------\n')
            output.write('tleap'.center(100)+'\n')
            output.write('----------------------------------------------------------------------------------------------------\n')
            output.write(result.stdout)
        #else, the forcefield is charmm and converting psf to mol2
        else:
            psf_to_mol2(keywords['pdb_file'])
            shutil.copy(keywords['pdb_file'], 'prot_autocap_fixed.pdb')
        
        #checking for integer charge in the mol2
        resi_output = check_resi_charges('prot_autocap_fixed.mol2')
        if resi_output[0] == 0:
            print(resi_output[1])
            sys.exit()
        #combining information from the cx pdb, protein pdb, and mol2 into a dataframe for easy handling
        #updating water charges if needed
        if 'ep_charge' in keywords:
            combine_data(keywords['o_charge'], keywords['h_charge'], keywords['ep_charge'])
        elif 'h_charge' in keywords:
            combine_data(keywords['o_charge'], keywords['h_charge'])
        else:
            combine_data()

        #checking created dataframe for integer charges
        df_output = check_df_charges()
        if df_output[0] ==0:
            print(df_output[1])
            sys.exit() 
        output.close()

        
        #if the cutoff is zero, creating dictionary with all atoms in the MM region
        if keywords['cutoff'] =='0':
            dictionary_nocut()
            #shutil.move('dictionary.dat', new_dir)
            shutil.copy(keywords['pdb_file'], f'{new_dir}/external.pdb')
        #elif the a template path has been specified, mapping the QM region to the QM region of the template
        elif 'template_path' in keywords:
            convert_dictionary(keywords['cutoff'], keywords['template_path'])
        #if there is no template specified then cutting the protein
        else:
            run_cut_protein('cx_autocap_fixed.pdb', seed, keywords['cutoff'])
            os.mkdir(f'data')
            shutil.move('QM.pdb', f'data/QM-sub-cut-protein-fragment-ligand.pdb')
            #if the protein was cut then move M3 atoms that are in different residues than the M2 atoms into the MM region
            move_m3s()
            #shutil.move('external.pdb', new_dir)
            #shutil.move('pre-dictionary.dat', new_dir)
            #shutil.move('ligand.pdb', new_dir)
        
        #if any cuts were made, cap the cut QM bonds with link hydrogens
        if keywords['cutoff']!='0':
            if 'amber_ff' in keywords:
                run_cap(ff_type = 'amber', path_to_env = keywords['env_path'])
            else:
                run_cap(ff_type = 'charmm', rtf = keywords['charmm_rtf'], prm = keywords['charmm_prm'])
            #redistribute charge based on charge scheme and write QM input file
            write_input(input_file, f'{new_dir}_psi4_file.py')
            if 'do_fsapt' in keywords:
                if keywords['do_fsapt'] == 'false':
                    write_QM(keywords['charge_scheme'], keywords['ligand_charge'], keywords['basis_set'], keywords['method'], f'{new_dir}_psi4_file.py', False)
            else:
                write_QM(keywords['charge_scheme'], keywords['ligand_charge'], keywords['basis_set'], keywords['method'], f'{new_dir}_psi4_file.py')

            #check the charges and number of atoms in the written QM input file
            check_QM_file()
        #write fsapt files
        if keywords['fisapt_partition'] == 'true':
            partition('CAPPED_qm.pdb')
        print(f"\u2728Sparcle-QC has sparkled\u2728")
    except Exception as e:
        print(f'\n An error has occured: {e}')
    finally:
        stop_flashing.set()
    
    # Wait for the flashing thread to finish
        flashing_thread.join()
