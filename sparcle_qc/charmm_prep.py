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

def fix_numbers_charmm(pdb_file: str) -> None:
    """
    When given a pdb, creates a new copy {pdb_file}_fixed.pdb that has the protein residues followed by waters and then the ligand
    Corrects for any mistakes in atom or residue numbering that may have been caused by manipulation of the system in pymol
    Ensures that the ligand atoms are labeled as HETATM

    Parameters
    ----------
    pdb_file: str
        path to pdb

    Returns
    -------
    None
    """

    charmmsys = pmd.load_file(pdb_file)
    with open('ligand.pdb') as lig:
        lig_lines = lig.readlines()
    lig_name = lig_lines[3][16:20].strip()
    
    out = open(f'{pdb_file[:-4]}_fixed.pdb', 'w') 
    with open(pdb_file) as w:
        lines = w.readlines()
    resnum =0
    atomnum = 0
    ligand_lines = []
    HOH_lines = []
    oldres = ''
    for line in lines:
        if 'HOH' not in line and 'TIP' not in line and len(line)>70 and line[16:20].strip() !=lig_name and (line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM'):
            atomnum +=1
            if line[22:26].strip()!=oldres:
                resnum+=1
                oldres = line[22:26].strip()
            atomtype = charmmsys.atoms[int(line[6:11].strip())-1].element_name
            out.write(f'ATOM  {atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{atomtype:>12}\n')
        elif 'HOH' in line or 'TIP' in line and (line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM'):
            HOH_lines.append(line)
        elif lig_name in line and (line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM'):
            ligand_lines.append(line)
        else:
            pass
    
    for line in HOH_lines:
        if len(line)>70:
            atomnum +=1
            if line[22:26].strip()!=oldres:
                resnum+=1
                oldres = line[22:26].strip()
            atomtype = charmmsys.atoms[int(line[6:11].strip())-1].element_name
            out.write(f'ATOM  {atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{atomtype:>12}\n')
    for line in ligand_lines:
        if len(line)>70 and line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM':
            atomnum +=1
            if line[22:26].strip()!=oldres:
                resnum+=1
                oldres = line[22:26].strip()
            atomtype = charmmsys.atoms[int(line[6:11].strip())-1].element_name
            out.write(f'HETATM{atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{atomtype:>12}\n')
    if 'cx' in pdb_file:
        out.write('CONECT\n')
    out.write('END')
    out.close()

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
    When given a CHARMM psf, converts the information encoded into the style of a mol2 

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

#this is not a charmm specific function, but is located in this module because of its dependence on parmed
def dictionary_nocut(cx_pdb:str = 'cx_autocap_fixed.pdb') -> None:
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
