import sys
import os
#parmed has a warning that doesn't need to be displayed to the terminal
with open(os.devnull, 'w') as devnull:
    old_stdout = sys.stdout
    sys.stdout = devnull
    try:
        import parmed as pmd
    finally:
        sys.stdout = old_stdout
import json

def combine_charmm(prot_file: str) -> None:
    """
    When given a CHARMM protein pdb, uses parmed to combine it with the ligand into a single complex pdb

    Parameters
    ----------
    pdb_file: str
        path to protein pdb

    Returns
    -------
    None
    """
    
    charmmprot = pmd.load_file(prot_file)
    charmmlig = pmd.load_file('ligand.pdb')
    structure = charmmprot+charmmlig
    structure.save('cx_autocap.pdb')

def psf_to_mol2(original_pdb: str) -> None:
    """
    When given a CHARMM psf, converts the information encoded into the style of an AMBER mol2 

    Parameters
    ----------
    pdb_file: str
        path to original pdb from the input file

    Returns
    -------
    None
    """
    mol2_path = 'prot_autocap_fixed.mol2'
    psf_path = original_pdb.replace('pdb', 'psf')
    pdb_path = 'cx_autocap_fixed.pdb'  
    coord_dict = {}
    ext = False
    with open(pdb_path) as pdb_file:
        pdb_lines = pdb_file.readlines()
    for line in pdb_lines:
        if line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM':
            coord_dict[line[6:11].strip()] = [line[30:38].strip(),line[38:46].strip(),line[46:54].strip()]
    
    with open(psf_path) as psf_file:
        psf_lines = psf_file.readlines()
    for num, line in enumerate(psf_lines):
        if '!NATOM' in line:
            num_atom = num
    
        if '!NBOND' in line:
            num_bond = num
            break
    atom_info = psf_lines[num_atom+1:num_bond-1]
    
    mol2_file = open(mol2_path, 'w')
    
    mol2_file.write('@<TRIPOS>MOLECULE\n')
    mol2_file.write(f'default_name\n ****  ****  ****  ****  ****\n')
    
    mol2_file.write('SMALL\n')
    mol2_file.write('USER_CHARGES\n')
    
    mol2_file.write('@<TRIPOS>ATOM\n')
    
    for line in atom_info:
        atom_id = line.split()[0]
        atom_name = line.split()[4]
        atom_type = line.split()[5]
        subst_id = line.split()[2]
        subst_name = line.split()[3]
        charge = line.split()[6]
        x = float(coord_dict[atom_id][0])
        y = float(coord_dict[atom_id][1])
        z = float(coord_dict[atom_id][2])
        charge = float(charge)
        status_bit ='****'
    
        mol2_file.write(f'{atom_id} {atom_name:<4} {x:>12.6f} {y:>12.6f} {z:>12.6f} {atom_type:>4} {subst_id:>6} {subst_name:>4} {charge:>11.4f} {status_bit}\n')
    
    mol2_file.write('@<TRIPOS>BOND\n')
    mol2_file.close()


def dictionary_nocut(cx_pdb = 'cx_autocap_fixed.pdb': str) -> None:
    """
    if the cutoff specified in the input file is 0 angstroms, the the entirety of the protein should be located in the MM region
    this function uses parmed to loop through each atom in the complex pdb and add into the MM list of atoms in dictionary.dat

    Parameters
    ----------
    pdb_file: str
        path to complex pdb

    Returns
    -------
    None
    """
    lig = pmd.load_file('ligand.pdb')
    lig_name = lig.residues[0].name
    d = {'MM': []}
    charmmprot = pmd.load_file(cx_pdb)
    for atom in charmmprot.atoms:
        if atom.residue.name != lig_name:
            d['MM'].append(atom.idx+1)
    with open('dictionary.dat', 'w+') as wfile:
        json.dump(d, wfile)
