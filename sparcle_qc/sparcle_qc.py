import sys
import traceback
import shutil
import subprocess
import os
import threading
import itertools
import time
import ast
from typing import Dict

from sparcle_qc.amber_prep import write_cpptraj, write_cpptraj_skip_autocap, write_tleap, autocap, skip_autocap, reorder_atoms_amber
from sparcle_qc.charmm_prep import psf_to_mol2, get_cx_pdb, reorder_atoms_charmm
from sparcle_qc.complex_tools import check_df_charges, check_mol2_charges, convert_seed, closest_contact
from sparcle_qc.combine_data import create_csv
from sparcle_qc.cut_protein import run_cut_protein 
from sparcle_qc.convert_dict import convert_dictionary 
from sparcle_qc.move_M3s import move_m3s
from sparcle_qc.create_est_inp import make_monomers, check_est_file, copy_input, write_est_file, ghost
from sparcle_qc.cap import run_cap
from sparcle_qc.make_fsapt_partition import fsapt_partition 
from sparcle_qc.make_frag_dirs import make_dirs

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
                split_line = line.split(':', 1)
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
                        if cutoff <= 0:
                            raise ValueError("Error: Cutoff radius must be greater than zero")
                    except ValueError as e:
                        print(e)
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
                if key_word == 'do_fsapt':
                    value = value.lower()
                    if value != 'true' and value != 'false':
                        print("Error: Invalid input file. do_fsapt is not true or false")
                        sys.exit()
                if key_word == 'amber_ff':
                    #potentially check if the ff is in amber
                    pass
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
                if key_word == 'software':
                    value = value.lower()
                    if value not in ['psi4','nwchem', 'q-chem']:
                        print("Error: Software not supported. Choose psi4, nwchem, or q-chem.")
                        sys.exit()
                if key_word == 'cp':
                    value = value.lower()
                    if value != 'true' and value != 'false':
                        print("Error: Invalid input file. cp is not true or false")
                        sys.exit()
                keywords[key_word]=value
    if 'software' not in keywords.keys():
        print('Error: Software not specified. Choose psi4, nwchem, or q-chem.')
        sys.exit()
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
        print('Error: Invalid input file. No charge scheme is provided')
        sys.exit()
    if 'method' not in keywords.keys():
        print('Error: Invalid input file. No method is provided')
        sys.exit()
    if 'basis_set' not in keywords.keys():
        print('Error: Invalid input file. No basis set is provided')
        sys.exit()
    if 'amber_ff' not in keywords.keys() and 'charmm_rtf' not in keywords.keys():
        print('Error: Invalid input file. CHARMM or Amber FF is not provided')
        sys.exit()
    if 'amber_ff' in keywords.keys():
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
    if 'fisapt_partition' not in keywords.keys():
        keywords['fisapt_partition'] = 'false'
    if 'do_fsapt' not in keywords.keys():
        keywords['do_fsapt'] = None
    if 'ep_charge' in keywords.keys():
        try:
            o_charge = keywords['o_charge']
            h_charge = keywords['h_charge']
        except KeyError:
            print('Error: Invalid input file. Oxygen and/or Hydrogen charges are not provided for water')
            sys.exit()
    if 'mem' not in keywords.keys():
        keywords['mem'] = '32 GB'
    if 'nthreads' not in keywords.keys():
        keywords['nthreads'] = '1'
    if 'cp' not in keywords.keys():
        keywords['cp'] = 'true'
    if 'nwchem_scratch' not in keywords.keys() and keywords['software'].lower() == 'nwchem':
        print('Error: nwchem_scratch not provided.')
        sys.exit()
    if 'nwchem_scratch' not in keywords.keys():
        keywords['nwchem_scratch'] = None
    if 'nwchem_perm' not in keywords.keys() and keywords['software'].lower() == 'nwchem':
        print('Error: nwchem_perm not provided')
        sys.exit()
    if 'nwchem_perm' not in keywords.keys():
        keywords['nwchem_perm'] = None
    if 'nwchem_scf' in keywords.keys():
        keywords['nwchem_scf'] = ast.literal_eval(keywords['nwchem_scf'])
        if isinstance(keywords['nwchem_scf'], dict) is False:
            print('Error: nwchem_scf is not a dictionary')
            sys.exit()
    else:
        keywords['nwchem_scf'] = None
    if 'nwchem_dft' in keywords.keys():
        keywords['nwchem_dft'] = ast.literal_eval(keywords['nwchem_dft'])
        if isinstance(keywords['nwchem_dft'], dict) is False:
            print('Error: nwchem_dft is not a dictionary')
            sys.exit()
    else:
        if keywords['method'].lower() == 'dft':
            keywords['nwchem_dft'] = {'xc':'b3lyp'}
        else:
            keywords['nwchem_dft'] = None
    if 'psi4_options' in keywords.keys():
        keywords['psi4_options'] = ast.literal_eval(keywords['psi4_options'])
        if isinstance (keywords['psi4_options'], dict) is False:
            print('Error: psi4_options is not a dictionary')
            sys.exit()
    else:
        keywords['psi4_options'] = {}
    if 'freeze_core' not in (key.lower() for key in keywords['psi4_options'].keys()):
        keywords['psi4_options']['freeze_core'] = 'true'
    if 'scf_type' not in (key.lower() for key in keywords['psi4_options'].keys()):
        keywords['psi4_options']['scf_type'] = 'df'
    if 'qchem_options' in keywords.keys():
        keywords['qchem_options'] = ast.literal_eval(keywords['qchem_options'])
        if isinstance (keywords['qchem_options'], dict) is False:
            print('Error: qchem_options is not a dictionary')
            sys.exit()
    elif keywords['software'].lower() == 'q-chem':
        keywords['qchem_options'] = {}
        if 'jobtype' not in (key.lower() for key in keywords['qchem_options'].keys()):
            if 'sapt' in keywords['method']:
                keywords['qchem_options']['JOBTYPE'] = 'xsapt'
            else:
                keywords['qchem_options']['JOBTYPE'] = 'sp'
    else:
        keywords['qchem_options'] = None
    if 'qchem_sapt' in keywords.keys():
        keywords['qchem_sapt'] = ast.literal_eval(keywords['qchem_sapt'])
        if isinstance (keywords['qchem_sapt'], dict) is False:
            print('Error: qchem_sapt is not a dictionary')
            sys.exit()
    else:
        keywords['qchem_sapt'] = {}
    if keywords['method'].lower() == 'sapt0' and keywords['software'] == 'q-chem':
        if 'algorithm' not in (key.lower() for key in keywords['qchem_sapt'].keys()):
            keywords['qchem_sapt']['algorithm'] = 'ri-mo'
        if 'basis' not in (key.lower() for key in keywords['qchem_sapt'].keys()):
            keywords['qchem_sapt']['basis'] = 'dimer'
    else:
        keywords['qchem_sapt'] = None
    if keywords['software'] == 'nwchem' and 'sapt' in keywords['method']:
        print('Error: SAPT is not available in NWChem. Choose a different method.')
        sys.exit()
    if 'other_amber_ff' in keywords.keys():
        try:
            ast.literal_eval(keywords['other_amber_ff'])
        except:
            print("Error: other_amber_ff is not a list of strings")
            sys.exit()
        keywords['other_amber_ff'] = ast.literal_eval(keywords['other_amber_ff'])
    else:
        keywords['other_amber_ff'] = []
    print(f"\u2728Sparcle-QC is sparkling\u2728\nBeginning file preparation for an embedded QM calculation of {keywords['pdb_file']} ")
    
    return keywords


def run_sparcle(input_file= None, user_options = None):
    """ 
    Given an input file, parses the specified options into a
    dictionary of keywords (or just takes in a dictionary) 
    and runs the necessary sparcle_qc functions to
    create a quantum chemistry software input file. Steps include
    preparing PDBs, obtaining MOL2s, carving out the QM region, capping
    the cut bonds with hydrogens, and writing an input file.

    Parameters
    ----------
    input_file: str
        path to input file
    user_options: Dict
        User provided dictionary with sparcle_qc parameters instead of providing an input file

    Returns
    -------
    number of QM atoms, number of MM atoms, charge of QM region, charge of MM region: List or Dict
        Number of atoms and charge of each region in a List for a SAPT computation or a Dictionary for a supermolecular computation with keys 'Complex', 'Protein', and 'Ligand'
        
    """
    
    #starting sparkle on command line
    flashing_thread = threading.Thread(target=flashing_sparkle)
    flashing_thread.start()
    try:
        #parsing input file into dictionary
        if input_file == None and user_options == None:
            print("Error: Input file or Dictionary not provided")
            sys.exit()
        if input_file != None and user_options != None:
            print("Error: Input file and Dictionary provided. Choose one.")
            sys.exit()
        if user_options != None:
            if 'input_filename' not in user_options: 
                print("Error: Input file name not provided")
                sys.exit()
            if '.' not in user_options['input_filename']: 
                print("Error: Input file extension not provided")
                sys.exit()
            try:
                str(user_options['input_filename'])
            except:
                print("Error: Input file not a string")
                sys.exit()
            with open(user_options['input_filename'], 'w') as inp:
                for key in user_options:
                    if key != 'input_filename':
                        inp.write(f'{key}: {user_options[key]}\n')
            input_file = user_options['input_filename']

        keywords = input_parser(input_file)
        #creating new directory for the created files, changing working directories, and copying necessary files into the new directory
        new_dir = '.'.join(input_file.split('.')[:-1])
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
            reorder_atoms_amber('prot_autocap.pdb')
            reorder_atoms_amber('cx_autocap.pdb')
            os.remove('prot_autocap.pdb')
        #otherwise, the forcefield is charmm and combining protein and ligand into a complex pdb
        else:
            get_cx_pdb(keywords['pdb_file'])
            reorder_atoms_charmm('cx_autocap.pdb')
        os.remove('cx_autocap.pdb')

        #obtaining seed containing information of which group to grow the QM region from
        if keywords['seed'] =='ligand':
            seed = 'ligand'
            seed_coords = 'ligand.pdb'
        else:
            seed, seed_coords = convert_seed(keywords['seed'], keywords['seed_file'])
        #if forcefield is amber, writing and running tleap to create mol2
        if 'amber_ff' in keywords:
            write_tleap(keywords['amber_ff'], keywords['water_model'], keywords['other_amber_ff'])
            result = subprocess.run(['tleap -f tleap.in'], text = True, shell = True, capture_output = True)
            output.write('----------------------------------------------------------------------------------------------------\n')
            output.write('tleap'.center(100)+'\n')
            output.write('----------------------------------------------------------------------------------------------------\n')
            output.write(result.stdout)
        #else, the forcefield is charmm and converting psf to mol2
        else:
            psf_to_mol2(keywords['pdb_file'])
            shutil.copy(keywords['pdb_file'], 'prot_autocap_fixed.pdb')
        min_dist = closest_contact('prot_autocap_fixed.pdb', seed_coords)
        if float(keywords['cutoff']) < min_dist:
            print(f'Error: Cutoff is less than the shortest intermolecular distance between the seed atom(s) and the protein. Please choose a value greater than {min_dist:.2f} Ang.')
            sys.exit()
        
        #checking for residue integer charges in the mol2
        resi_output = check_mol2_charges('prot_autocap_fixed.mol2')
        if resi_output[0] == 0:
            print(resi_output[1])
            sys.exit()
        #combining information from the cx pdb, protein pdb, and mol2 into a dataframe for easy handling
        #updating water charges if specified by user
        if 'ep_charge' in keywords:
            create_csv(keywords['o_charge'], keywords['h_charge'], keywords['ep_charge'])
        elif 'h_charge' in keywords:
            create_csv(keywords['o_charge'], keywords['h_charge'])
        else:
            create_csv()

        #checking created dataframe for residue integer charges
        df_output = check_df_charges()
        if df_output[0] ==0:
            print(df_output[1])
            sys.exit() 
        output.close()

        
        #if the a template path has been specified, mapping the QM region from the QM region of the template
        if 'template_path' in keywords:
            convert_dictionary(keywords['cutoff'], keywords['template_path'])
        # if resis_per_fragment is specified, cut protein into many fragments with N residues
        elif 'resis_per_fragment' in keywords:
            run_cut_protein('cx_autocap_fixed.pdb', seed, keywords['cutoff'], keywords['resis_per_fragment'])
            dirs_list = make_dirs()
            for frag_dir in dirs_list:
                os.chdir(frag_dir)
                move_m3s()
                #cap the cut QM bonds with link hydrogens
                if 'amber_ff' in keywords:
                    run_cap(ff_type = 'amber')
                else:
                    run_cap(ff_type = 'charmm', rtf = keywords['charmm_rtf'], prm = keywords['charmm_prm'])
                os.chdir('../')
            print(dirs_list)
            exit()
        #if there is no template specified or resis_per_fragment then cutting the protein
        else:
            run_cut_protein('cx_autocap_fixed.pdb', seed, keywords['cutoff'])
            os.mkdir(f'data')
            shutil.move('QM.pdb', f'data/QM-sub-cut-protein-fragment-ligand.pdb')
            #for each cut, ensure frontier MM atoms are part of one unique residue
            move_m3s()
            #shutil.move('external.pdb', new_dir)
            #shutil.move('pre-dictionary.dat', new_dir)
            #shutil.move('ligand.pdb', new_dir)
        
        #cap the cut QM bonds with link hydrogens
        if 'amber_ff' in keywords:
            run_cap(ff_type = 'amber')
        else:
            run_cap(ff_type = 'charmm', rtf = keywords['charmm_rtf'], prm = keywords['charmm_prm'])
        #redistribute charge based on charge scheme and write QM input file
        qm_lig, c_QM, qm_pro, mm_env = make_monomers(keywords['charge_scheme'])
        ext = {'psi4':'.py', 'nwchem':'.in', 'q-chem':'.in'}
        sapt_inp_filename = f'{new_dir}_' + keywords['software'] + '_file' + ext[keywords['software']]
        if 'sapt' in keywords['method'].lower():
            copy_input(input_file, sapt_inp_filename, keywords['software'])
            write_est_file(keywords['software'], qm_lig, c_QM, qm_pro, '', mm_env, sapt_inp_filename, keywords['ligand_charge'], keywords['method'], keywords['basis_set'], keywords['mem'], keywords['nthreads'], keywords['do_fsapt'], keywords['nwchem_scratch'], keywords['nwchem_perm'], keywords['nwchem_scf'], keywords['nwchem_dft'], keywords['psi4_options'], keywords['qchem_options'], keywords['qchem_sapt'])
#            #check the charges and number of atoms in the written QM input file
            qm_atoms, mm_atoms, qm_charge, mm_charge = check_est_file(sapt_inp_filename)
        else:
            if keywords['cp'] == 'true':
                ghost_lig, lig_uniq_elements = ghost(qm_lig, keywords['software'])
                ghost_pro, prot_uniq_elements = ghost(qm_pro, keywords['software'])
                ghost_charge = 0
            else:
                ghost_lig = None
                ghost_pro = None
                lig_uniq_elements = None
                prot_uniq_elements = None
                ghost_charge = None
            lig_inp_filename = f'{new_dir}_' + keywords['software'] + '_file_lig' + ext[keywords['software']]
            copy_input(input_file, lig_inp_filename, keywords['software'])
            write_est_file(keywords['software'], qm_lig, ghost_charge, ghost_pro, prot_uniq_elements, None, lig_inp_filename, keywords['ligand_charge'], keywords['method'], keywords['basis_set'], keywords['mem'], keywords['nthreads'], False, keywords['nwchem_scratch'], keywords['nwchem_perm'], keywords['nwchem_scf'], keywords['nwchem_dft'], keywords['psi4_options'], keywords['qchem_options'])
            prot_inp_filename = f'{new_dir}_' + keywords['software'] + '_file_prot' + ext[keywords['software']]
            copy_input(input_file, prot_inp_filename, keywords['software'])
            write_est_file(keywords['software'], ghost_lig, c_QM, qm_pro, lig_uniq_elements, mm_env, prot_inp_filename, ghost_charge, keywords['method'], keywords['basis_set'], keywords['mem'], keywords['nthreads'], None, keywords['nwchem_scratch'], keywords['nwchem_perm'], keywords['nwchem_scf'], keywords['nwchem_dft'], keywords['psi4_options'], keywords['qchem_options'])
            cx_inp_filename = f'{new_dir}_' + keywords['software'] + '_file_cx' + ext[keywords['software']]
            copy_input(input_file, cx_inp_filename, keywords['software'])
            write_est_file(keywords['software'], qm_lig, c_QM, qm_pro, None, mm_env, cx_inp_filename, keywords['ligand_charge'], keywords['method'], keywords['basis_set'], keywords['mem'], keywords['nthreads'], None, keywords['nwchem_scratch'], keywords['nwchem_perm'], keywords['nwchem_scf'], keywords['nwchem_dft'], keywords['psi4_options'], keywords['qchem_options'])
            qm_atoms_lig, mm_atoms_lig, qm_charge_lig, mm_charge_lig = check_est_file(lig_inp_filename)
            qm_atoms_pro, mm_atoms_pro, qm_charge_pro, mm_charge_pro = check_est_file(prot_inp_filename)
            qm_atoms_cx, mm_atoms_cx, qm_charge_cx, mm_charge_cx = check_est_file(cx_inp_filename)

        #write fsapt files
        if keywords['fisapt_partition'] == 'true':
            fsapt_partition('CAPPED_qm.pdb')
        
        print(f"\u2728Sparcle-QC has sparkled\u2728")
        if 'sapt' in keywords['method'].lower():
            return qm_atoms, mm_atoms, qm_charge, mm_charge
        else:
            return {'Complex':[qm_atoms_cx, mm_atoms_cx, qm_charge_cx, mm_charge_cx], 'Ligand':[qm_atoms_lig, mm_atoms_lig, qm_charge_lig, mm_charge_lig], 'Protein':[qm_atoms_pro, mm_atoms_pro, qm_charge_pro, mm_charge_pro]}
    except Exception as e:
        error_type = type(e).__name__
        error_message = str(e)
        error_traceback = traceback.format_exc()

        print(f'\nAn error has occurred:')
        print(error_traceback)
    finally:
        cur_path = os.getcwd()
        relative_dir = cur_path.split('/')[-1]
        try:
            if relative_dir==new_dir:
                os.chdir('../')
        except:
            pass
        stop_flashing.set()
    
    # Wait for the flashing thread to finish
        flashing_thread.join()
    
def main():
    if len(sys.argv) != 2:
        print("Usage: sparcle_qc <input_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    # Your main script logic here
    run_sparcle(input_file)

if __name__ == "__main__":
    main()
