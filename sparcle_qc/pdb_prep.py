import numpy as np
import pandas as pd
from typing import Tuple

def fix_numbers(pdb_file: str) -> None:
    """
    When given a pdb, creates a new copy {pdb_file}_fixed.pdb that has the protein residues followed by waters and then the ligand
    Corrects for any mistakes in atom or residue numbering

    Parameters
    ----------
    pdb_file: str
        path to pdb

    Returns
    -------
    None
    """

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
            out.write(f'ATOM  {atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{line[11:16].strip()[0]:>12}\n')
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
            out.write(f'ATOM  {atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{line[66:78].strip():>12}\n')
    for line in ligand_lines:
        if len(line)>70 and line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM':
            atomnum +=1
            if line[22:26].strip()!=oldres:
                resnum+=1
                oldres = line[22:26].strip()
            out.write(f'HETATM{atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{line[66:78].strip():>12}\n')
    if 'cx' in pdb_file:
        out.write('CONECT\n')
    out.write('END')
    out.close()

def get_coords(pdb_file1: str, atom_id1: str) -> Tuple[str, str, str]:
    """
    extracts the x,y,z coordinates of a specified atom in a specified pdb

    Parameters
    ----------
    pdb_file: str
        path to uncapped complex pdb
    
    Returns
    -------
    x_coord, y_coord, z_coord: str, str, str
        coordinates of specified atom
    """

    with open(pdb_file1, 'r') as file:
        lines = file.readlines()
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
                if l[6:11].strip() == atom_id1:
                    x_coord = l[29:38].strip()
                    y_coord = l[38:46].strip()
                    z_coord = l[46:54].strip()
                    return x_coord, y_coord, z_coord

def match_coords(x:str, y:str, z:str, pdb_file2:str) -> str:
    """
    returns the atom number that corresponds to the specified x,y,z coordinates in the specified pdb_file

    Parameters
    ----------
    pdb_file: str
        path to uncapped complex pdb

    Returns
    -------
    atom_id: str
        atom id that corresponds to the specified coordinates in the specified pdb file
    """
    with open(pdb_file2, 'r') as file:
        lines = file.readlines()
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
               x_coord = l[29:38].strip()
               y_coord = l[38:46].strip()
               z_coord = l[46:54].strip()
               if x_coord == x and y_coord == y and z_coord == z:
                   return l[6:11].strip()
def convert_atom_id(seed: str, seed_file:str, new_pdb ='cx_autocap_fixed.pdb':str) -> str:
    """
    finds the atom the corresponds to atom number = seed in seed_file pdb and then maps this to the newly created complex pdb

    Parameters
    ----------
    seed: str
        atom number in seed_file
    seed_file: str
        file to original pdb that contains seed
    new_pdb: str
        path to new pdb

    Returns
    -------
    cx_id: str
        atom in new_pdb that corresponds to the atom seed in seed_file
    """
    x1,y1,z1 = get_coords(seed_file, seed)
    cx_id = match_coords(x1, y1, z1, new_pdb)
    return cx_id

def check_resi_charges(mol2_file: str) -> Tuple[int, str]:
    """
    checks the sum of the point charges for each residue listed in the mol2 and ensures that it sums to an integer

    Parameters
    ----------
    mol2_file: str
        path to mol2 file

    Returns
    -------
    error_message: tuple(int, str)
        tuple where the first entry is an integer. 0 = at least one residue is not integer charge. 1 = all residues are integer charge.
        the second entry is a string containing information about which residue is not integer charge and the corresponding fractional charge
    """
    total_charge = 0
    tolerance = 0.0001
    failed = False
    with open(mol2_file, 'r') as file:
        all_lines = file.readlines()
        for n,line in enumerate(all_lines):
            if 'ATOM' in line:
                start = n
            elif 'BOND' in line:
                end = n
        lines = all_lines[start+1:end]
        init_resi = lines[0].split()[-3]
        charge = 0
        for line in lines:
            resi_name = line.split()[-3]
            resi_number = line.split()[-4]
            if resi_number == init_resi:
                charge += float(line.split()[-2])
            else:
                total_charge += charge
                if np.abs(np.round(charge, 0) - charge) > tolerance:
                    return (0, f'Residue {init_resi} charge: {charge}. Fix pdb/mol2 such that this residue has integer charge and restart')
                    failed = True
                init_resi = resi_number
                charge = float(line.split()[-2])
    if not failed:
        return (1, 'passed')

def check_df_charges() -> Tuple[int, str]:
    """
    checks the sum of the point charges for each residue listed in the dataframe and ensures that it sums to an integer

    Parameters
    ----------

    Returns
    -------
    error_message: tuple(int, str)
        tuple where the first entry is an integer. 0 = at least one residue is not integer charge. 1 = all residues are integer charge.
        the second entry is a string containing information about which residue is not integer charge and the corresponding fractional charge
    """
    total_charge = 0
    tolerance = 0.0001 
    failed = False
    
    charges = {}
    with open('dataframe.csv', 'r') as file:
        lines = file.readlines()[1:]
        charge = 0
        for line in lines:
            resi_name = line.split(',')[2]
            if resi_name in charges.keys():
                charges[resi_name] += float(line.split(',')[-2])
            else:
                charges[resi_name] = float(line.split(',')[-2])
    
    total_charge = 0
            
    df = pd.read_csv('dataframe.csv', index_col='PDB_ID')
    
    # loop through residues and their sum of point charges
    for k,v in charges.items():
        total_charge += v
        # check for integer charge
        if np.abs(np.round(v, 0) - v) > tolerance:
            failed = True
            return (0, f'Residue {k} charge: {v}. Fix pdb/mol2 such that this residue has integer charge and restart')
        if 'WAT' in k or 'HOH' in k:
            # check for neutral waters
            if v != 0:
                failed = True
                return (0, f'Water charge ({k}) does not equal zero. Please choose a new water model or add custom charges in the input. See dataframe.csv for current charges.')
        if np.round(v,0) != 0:
            # for non-zero residue charges:
            int_q = np.round(v, 0)
            if int_q < 0:
                q_str = f'{int(np.abs(int_q))}' + '-'
            elif int_q > 0:
                q_str = f'{int(np.abs(int_q))}' + '+'
            idxs = df[df['PDB_RES'] == k].index.tolist()
            for idx in idxs:
                # check for formal charge already noted in df (1+, 1-, etc.)
                if '+' in df.at[idx, 'AT_LABEL'] or '-' in df.at[idx, 'AT_LABEL']:
                    df_charge = df.at[idx, 'AT_LABEL'][-2:]
                    if df_charge != q_str:
                        # make sure dataframe formal charge equals sum of point charges. 
                        return (1,f'PDB formal charge ({df_charge}) does not equal sum of partial charges in mol2 ({q_str}) for {k}. Proceeding with the formal charge of residue set to the sum of partial charges.')
    return (1,None)
