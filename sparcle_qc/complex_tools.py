import numpy as np
import pandas as pd
from typing import Tuple
import parmed as pmd


def id_to_coords(pdb_file1: str, atom_id1: str) -> Tuple[str, str, str]:
    """
    extracts the x,y,z coordinates of a specified atom in a specified pdb

    Parameters
    ----------
    pdb_file: str
        path to pdb
    
    Returns
    -------
    x_coord, y_coord, z_coord: Tuple[str, str, str]
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

def coords_to_new_id(x:str, y:str, z:str, pdb_file2:str) -> str:
    """
    returns the atom number that corresponds to the specified x,y,z
    coordinates in the specified pdb_file

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

def convert_atom_id(seed: str, seed_file:str, new_pdb:str ='cx_autocap_fixed.pdb') -> str:
    """
    finds the atom the corresponds to atom number = seed in seed_file
    pdb and then maps this to the newly created complex pdb

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
    coords: list
        coordinates of seed atom
    """
    x1,y1,z1 = id_to_coords(seed_file, seed)
    cx_id = coords_to_new_id(x1, y1, z1, new_pdb)
    coords = [float(x1),float(y1),float(z1)]
    return cx_id, coords

def check_mol2_charges(mol2_file: str) -> Tuple[int, str]:
    """
    checks the sum of the point charges for each residue listed in the
    mol2 and ensures that it sums to an integer

    Parameters
    ----------
    mol2_file: str
        path to mol2 file

    Returns
    -------
    error_message: Tuple[int, str]
        tuple where the first entry is an integer. 
        0 = at least one residue is not integer charge. 
        1 = all residues are integer charge.
        The second entry is a string containing information about which 
        residue is not integer charge and the corresponding fractional charge.
    """
    total_charge = 0
    #tolerance is how far away a residue can be from an integer without throwing an error
    #this was determined because amber lists charge to 3 decimal places, 
    #so anything different to the fourth decimal place or further is insignificant
    tolerance = 0.0001
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
                init_resi = resi_number
                charge = float(line.split()[-2])
    return (1, 'passed')

def check_df_charges() -> Tuple[int, str]:
    """
    checks the sum of the point charges for each residue listed in the
    dataframe and ensures that it sums to an integer

    Parameters
    ----------
    None

    Returns
    -------
    return_message: Tuple[int, str]
        tuple where the first entry is an integer. 
        0 = at least one residue is not integer charge. 
        1 = all residues are integer charge.
        The second entry is a string containing information about which 
        residue is not integer charge and the corresponding fractional charge.
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
            return_message = (0, f'Residue {k} charge: {v}. Fix pdb/mol2 such that this residue has integer charge and restart')
        if 'WAT' in k or 'HOH' in k:
            # check for neutral waters
            if v != 0:
                failed = True
                return_message = (0, f'Water charge ({k}) does not equal zero. Please choose a new water model or add custom charges in the input. See dataframe.csv for current charges.')
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
                        print(f'PDB formal charge ({df_charge}) does not equal sum of partial charges in mol2 ({q_str}) for {k}. Proceeding with the formal charge of residue set to the sum of partial charges.')
                        # if not, remove the formal charge listed
                        df.at[idx, 'AT_LABEL'] = df.at[idx, 'AT_LABEL'][:-2]
                        # then add the sum of partial charges to the beta carbon atom label.
                        cb_idx = df[(df['PDB_RES']  == k) & (df['PDB_AT'] == 'CB')].index.tolist()
                        if len(cb_idx) > 0:
                            df.at[cb_idx[0], 'AT_LABEL'] = str(df.at[cb_idx[0], 'AT_LABEL']) + q_str
                        else:
                            # if there is no beta carbon (likely a metal ion)
                            possible_idx = df[(df['PDB_RES']  == k)].index.tolist()
                            df.at[possible_idx[0], 'AT_LABEL'] = str(df.at[possible_idx[0], 'AT_LABEL']) + q_str
                else:
                    # if there is no formal charge already listed in df, add it to the beta carbon
                    cb_idx = df[(df['PDB_RES']  == k) & (df['PDB_AT'] == 'CB')].index.tolist()
                    if len(cb_idx) > 0:
                        if '+' not in df.at[cb_idx[0],'AT_LABEL'] and '-' not in df.at[cb_idx[0],'AT_LABEL']:
                            df.at[cb_idx[0], 'AT_LABEL'] = str(df.at[cb_idx[0], 'AT_LABEL']) + q_str
                    else:
                        possible_idx = df[(df['PDB_RES']  == k)].index.tolist()
                        if '+' not in df.at[possible_idx[0],'AT_LABEL'] and '-' not in df.at[possible_idx[0],'AT_LABEL']:
                            df.at[possible_idx[0], 'AT_LABEL'] = str(df.at[possible_idx[0], 'AT_LABEL']) + q_str
    
            df.to_csv('dataframe.csv')
    if not failed:
        return_message = (1, 'worked')

    return return_message

def closest_contact(pdb1, pdb2):
    """
    Calculate the distance of the closest contact between two PDBs, or one PDB and one coordinate.
    
    Args:
    pdb1 (np.ndarray): A PDB representing the first set of XYZ coordinates.
    pdb2 (np.ndarray): A PDB representing the second set of XYZ coordinates or list of coords of one atom (seed).
    
    Returns:
    float: The minimum distance between any two points from pdb1 and pdb2.
    """
    pmdpdb1 = pmd.load_file(pdb1)
    coords1 = pmdpdb1.coordinates
    if 'pdb' in pdb2:
        pmdpdb2 = pmd.load_file(pdb2)
        coords2 = pmdpdb2.coordinates
    else:
        coords2 = pdb2

    # Ensure the input arrays are NumPy arrays
    coords1 = np.array(coords1)
    coords2 = np.array(coords2)

    # Compute the distance matrix
    distances = np.linalg.norm(coords1[:, np.newaxis] - coords2, axis=2)

    # Find the minimum distance
    min_distance = np.min(distances)

    return min_distance
