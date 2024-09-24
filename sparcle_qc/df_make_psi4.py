import numpy as np
import json
import os
import pandas as pd
from glob import glob
from typing import List, Dict, Tuple

def qm_ligand(LIGAND_PDB_PATH: str) -> List[str]:
    """
    This functions converts a pdb into a list of lines in the pdb

    Parameters
    ----------
    LIGAND_PDB_PATH: str
        path to the ligand pdb

    Returns
    -------
    lines: List[str]
        list of the lines of the pdb
    """
    #gets QM ligand lines from ligand PDB file
    lines = []
    with open(LIGAND_PDB_PATH, 'r') as pdbfile:
        pdblines = pdbfile.readlines()
        for l in pdblines[:-1]:
            if l[0:6].strip() != 'CRYST1':
                if l[0:6].strip() != 'CONECT':
                    lines.append(l)
    return lines

def pdb_to_xyz(PDB_lines: List[str]) -> List[str]:
    """
    This functions converts a list pdb lines into a list of xyzs with atom types

    Parameters
    ----------
    PDB_lines: str
        path to the ligand pdb

    Returns
    -------
    xyz: List[str]
        list of the atom types along with xyz coords
    """
    # translates PDB lines to XYZ lines
    xyz = []
    for l in PDB_lines:
        if '-' in l[66:78].strip() or '+' in l[66:78].strip():
            xyz.append(l[66:78].strip()[:-2])    
        else:
            xyz.append(l[66:78].strip())
        x_coord = l[29:38].strip()
        y_coord = l[38:46].strip()
        z_coord = l[46:54].strip()
        xyz.append(x_coord)
        xyz.append(y_coord)
        xyz.append(z_coord+'\n')
    return xyz

def atoms_to_pdb_lines(CAPPED_PDB_PATH:str, atoms:List[str]) -> List[str]:
    """
    This functions converts a set of atoms in a PDB into List of information from those specified atoms in the PDB

    Parameters
    ----------
    CAPPED_PDB_PATH: str
        path to the capped complex pdb
    atoms: List[str]
        a list of atoms in the CAPPED_PDB_PATH to create a pdb for

    Returns
    -------
    lines: List[str]
        list of the lines of the pdb
    """
    # gets indexed lines from PDB file
    pdb_lines = []
    with open(CAPPED_PDB_PATH, 'r') as pdbfile:
        pdblines = pdbfile.readlines()
        for i in atoms:
            for l in pdblines[:-1]:
                if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM': 
                    if int(l[6:11].strip()) == int(i):
                        pdb_lines.append(l)
    return pdb_lines

def write_pdb(CAPPED_PDB_PATH:str, atoms:List[str], CAPPED_QM_PATH:str, LIGAND_PATH:str) -> None:
    """
    This functions converts a set of atoms in a PDB into a new PDB

    Parameters
    ----------
    CAPPED_PDB_PATH: str
        path to the capped complex pdb
    atoms: List[str]
        a list of atoms in the CAPPED_PDB_PATH to create a pdb for
    CAPPED_QM_PATH: str
        path to write the new pdb

    Returns
    -------
    None
    """
    # gets indexed lines from PDB file
    atomnum = 0
    output_pdb = open(CAPPED_QM_PATH, 'w')
    with open(LIGAND_PATH, 'r') as pdbfile:
        pdblines = pdbfile.readlines()
        for line in pdblines[:-1]:
            if line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM': 
                atomnum+=1
                output_pdb.write(f'{line[0:6].strip():<6}{atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{line[22:26].strip():>4}{line[26:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{line[66:78].strip():>12}\n')
    with open(CAPPED_PDB_PATH, 'r') as pdbfile:
        pdblines = pdbfile.readlines()
        for i in atoms:
            for line in pdblines[:-1]:
                if line[0:6].strip() == 'ATOM' or line[0:6].strip() == 'HETATM': 
                    if int(line[6:11].strip()) == int(i):
                        atomnum+=1
                        output_pdb.write(f'{line[0:6].strip():<6}{atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{line[22:26].strip():>4}{line[26:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{line[66:78].strip():>12}\n')
    output_pdb.close()
    return 

def MM_xyz_to_charge_array(MM_lines_c: List[str], MM_lines_0: List[str], MOL2_PATH:str, MM_for_array: List[str]) -> List[str]:
    #TODO should we always just initialize MM_for_array to be empty at the start of this function instead of passing it as a parameter
    """
    Given an array of MM coordinates, this function returns the array of coordinates and charges needed for the QM input file 

    Parameters
    ----------
    MM_lines_c: List[str]

    MM_lines_0: list[str]

    MOL2_PATH: str
        path to mol2

    MM_for_array: List[str]
        List of charges of each MM atom along with their XYZ coordinates

    Returns
    -------
    MM_for_array: List[str]
        list of the charges along with xyz coords
    """
    # takes MM lines from PDB file, matches them to MOL2 file, gets charges, then builds array of MM atoms
    with open(MOL2_PATH, 'r', encoding="iso-8859-1") as mol2file:
        H2_WAT_atms = []
        lines = mol2file.readlines()
        for n,l in enumerate(lines):                                        # filter lines for atoms only
            if 'ATOM' in l:
                start = n
            elif 'HETATM' in l:
                start = n
            elif 'BOND' in l:
                end = n
        mollines = lines[start+1:end]
        for l in mollines:                                           # storing the next three lines for getting OPC dummy atom
            mol2idx = int(l.split()[0])                                          #- getting mol2 atom id / index
            mol2at = str(l.split()[1])                                           #- atom type
            mol2resi = str(l.split()[7])                                         #- residue
            mol2x = float("{:.3f}".format(float(l.split()[2])))             # matching lines in mol2 to MM XYZ coords bc
            mol2y = float("{:.3f}".format(float(l.split()[3])))             # indices btwn mol2 and pdb don't match
            mol2z = float("{:.3f}".format(float(l.split()[4])))
            if MM_lines_c != None:
                for i in range(int(len(MM_lines_c)/4)):             # looping through each [atom,x,y,z]
                    mmx = float(MM_lines_c[1+(4*i)])
                    mmy = float(MM_lines_c[2+(4*i)])
                    mmz = float(MM_lines_c[3+(4*i)])
                    if mmx == mol2x and mmy == mol2y and mmz == mol2z:
                        charge = l.split()[-2]
                        MM_for_array.append(str(charge))
                        MM_for_array.append(str(mmx))
                        MM_for_array.append(str(mmy))
                        MM_for_array.append(str(mmz)+'\n')
                        if mol2at == 'H2' and mol2resi == 'WAT':
                            H2_WAT_atms.append(mol2idx)
            if MM_lines_0 != None:
                for i in range(int(len(MM_lines_0)/4)):                           # looping through each [atom,x,y,z]
                    mmx = float(MM_lines_0[1+(4*i)])
                    mmy = float(MM_lines_0[2+(4*i)])
                    mmz = float(MM_lines_0[3+(4*i)])
                    if mmx == mol2x and mmy == mol2y and mmz == mol2z:
                        charge = l.split()[-2]
                        MM_for_array.append(str(0.000))                           # zeroing out the charge
                        MM_for_array.append(str(mmx))
                        MM_for_array.append(str(mmy))
                        MM_for_array.append(str(mmz)+'\n')
            if H2_WAT_atms != None:
                for i in H2_WAT_atms:
                    if mol2idx == i+1:
                        charge = l.split()[-2]
                        MM_for_array.append(str(charge))
                        MM_for_array.append(str(mol2x))
                        MM_for_array.append(str(mol2y))
                        MM_for_array.append(str(mol2z)+'\n')
    return MM_for_array
 
def SEE_atoms(num_bonds_broken:int, with_HL:Dict[str,List[int]]) -> List[str]:
    """
    Returns all M1, M2, and M3 atoms

    Parameters
    ----------
    num_bonds_broken: int
        number of bonds broken that need to be capped
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region

    Returns
    -------
    MM_for_array: List[str]
        list of the charges along with xyz coords
    """
    MM_atoms = []
    for bond in range(1,num_bonds_broken+1):
        if f'M1_{bond}' in with_HL and f'M2_{bond}' in with_HL and f'M3_{bond}' in with_HL:
            MM_atoms += (with_HL[f'M1_{bond}'] + with_HL[f'M2_{bond}'] + with_HL[f'M3_{bond}'])
        elif f'M1_{bond}' in with_HL and f'M2_{bond}' in with_HL and f'M3_{bond}' not in with_HL:
            MM_atoms += (with_HL[f'M1_{bond}'] + with_HL[f'M2_{bond}'])
        else:
            MM_atoms += (with_HL[f'M1_{bond}'])
    return MM_atoms

def SEE(MOL2_PATH:str, CAPPED_PDB_PATH: str, with_HL: Dict[str,List[int]], num_bonds_broken: int) -> List[str]:
    """
    Makes no changes to the MM region

    Parameters
    ----------
    num_bonds_broken: int
        number of bonds broken that need to be capped
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region

    Returns
    -------
    MM_for_array: List[str]
        list of the charges along with xyz coords
    """
    seeatoms = SEE_atoms(num_bonds_broken, with_HL) #get atoms from indexing
    seelines = atoms_to_pdb_lines(CAPPED_PDB_PATH, seeatoms) #get PDB lines from atoms list
    seexyz = pdb_to_xyz(seelines) # get XYZ lines from PDB lines
    mm_env = MM_xyz_to_charge_array(seexyz, None, MOL2_PATH, MM_for_array) # get charges from XYZ lines
    return mm_env

def Z1_atoms_zero(num_bonds_broken: int, with_HL: Dict[str,List[int]]) -> List[int]:
    """
    returns the M1 atoms that need to be zeroed

    Parameters
    ----------
    num_bonds_broken: int
        number of bonds broken that need to be capped
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region

    Returns
    -------
    MM_atoms: List[int]
        list of the MM atoms to zero
    """
    MM_atoms = []
    for bond in range(1,num_bonds_broken+1):
        MM_atoms += with_HL[f'M1_{bond}']
    return MM_atoms

def Z1_atoms_charge(num_bonds_broken:int, with_HL: Dict[str,List[int]]) -> List[str]:
    """
    returns the M2 and M3 atoms 

    Parameters
    ----------
    num_bonds_broken: int
        number of bonds broken that need to be capped
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region

    Returns
    -------
    MM_atoms: List[int]
        list of the MM atoms to keep
    """
    MM_atoms = []
    for bond in range(1,num_bonds_broken+1):
        if f'M2_{bond}' in with_HL.keys():
            MM_atoms += with_HL[f'M2_{bond}'] 
        if f'M3_{bond}' in with_HL.keys():
            MM_atoms += with_HL[f'M3_{bond}']
    return MM_atoms

def Z1(df:pd.DataFrame, with_HL:Dict[str,List[int]], num_bonds_broken: int) -> List[str]:
    """
    zeroes the M1 atoms

    Parameters
    ----------
    df: pd.DataFrame
        pandas dataframe containing all atoms in the system
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region
    num_bonds_broken: int
        number of bonds broken that need to be capped

    Returns
    -------
    MM_for_array: List[str]
        list of the charges along with xyz coords
    """
    MM_for_array = []
    z1atoms0 = Z1_atoms_zero(num_bonds_broken, with_HL)
    z1atomsc = Z1_atoms_charge(num_bonds_broken, with_HL) + with_HL['MM']
    ext_with_wat = z1atoms0 + z1atomsc
    for idx in ext_with_wat:
        if df.loc[idx, 'MOL2_RES'] == 'WAT' and df.loc[idx, 'MOL2_AT'] == 'OW':
            df_idx = df.loc[idx, 'MOL2_ID'] + 3.5
            ext_with_wat.append(df_idx)
    ext_wat_df = df.loc[df.index.isin(ext_with_wat)]
    for idx in ext_wat_df.index:
        if idx in z1atoms0:
            ext_wat_df.at[idx, 'q'] = 0
    MM_df = ext_wat_df[['q', 'X', 'Y', 'Z']]
    for x in MM_df.values.tolist():
        for n,y in enumerate(x):
            if n < 3:
                MM_for_array.append(str(y))
            else:
                MM_for_array.append(str(y)+'\n')
    mm_env = MM_for_array
    return mm_env

def DZn(k:int, df:pd.DataFrame, with_HL:Dict[str,List[int]], num_bonds_broken:int) -> List[str]:
    """
    zeros the ___ atoms and redistributes that charge into the ___ atoms

    Parameters
    ----------
    df: pd.DataFrame
        pandas dataframe containing all atoms in the system
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region
    num_bonds_broken: int
        number of bonds broken that need to be capped

    Returns
    -------
    MM_for_array: List[str]
        list of the charges along with xyz coords
    """
    output = open(glob('*.out')[0], 'a')
    MM_for_array = []
    MM_used = [] # keep track of MM atoms placed in psi4 file already
    for bond in range(1, num_bonds_broken+1):
        output.write(f'WORKING ON BOND {bond}\n')
        M1_atom = with_HL[f'M1_{bond}']
        M1atom_at = df.at[M1_atom[0], 'PDB_AT']
        M1atom_res = df.at[M1_atom[0], 'PDB_RES']
        MM_atoms = []
        zeroed = []
        M2_atoms = with_HL[f'M2_{bond}']
        residues = [df.at[x, 'PDB_RES'] for x in M2_atoms]
        QM_atoms = with_HL[f'Q1_{bond}'] + with_HL['QM']
        if k == 1:
            zeroed += with_HL[f'M1_{bond}']
            MM_atoms += Z1_atoms_charge(num_bonds_broken, with_HL) + with_HL['MM']
        elif k == 2:
            zeroed += with_HL[f'M1_{bond}'] + with_HL[f'M2_{bond}']
            MM_atoms += Z2_atoms_charge(num_bonds_broken, with_HL) + with_HL['MM']
        elif k == 3:
            zeroed += with_HL[f'M1_{bond}'] + with_HL[f'M2_{bond}'] + with_HL[f'M3_{bond}']
            MM_atoms += with_HL['MM']
        for x in zeroed:
            MM_used.append(x)
        # get MM atoms involved in redistribution
        MMdf = df.loc[df.index.isin(MM_atoms)]
        MMdf.to_csv('MMdfres.csv')
        dist_df = MMdf[MMdf['PDB_RES'].isin(residues)] # get all MM atoms in involved residues
        for x in dist_df.index:
            MM_used.append(x)
        # calculate charge to redistribute
        QMdf = df.loc[df.index.isin(QM_atoms)]
        QM_res_df = QMdf[QMdf['PDB_RES'].isin(residues)] 
        QM_res_charge = QM_res_df['q'].sum()
        group_df = df[df['PDB_RES'].isin([M1atom_res])] 
        group_charge = group_df['q'].sum()
        if M1atom_at == 'CA':
            q_to_dist = QM_res_charge
        elif M1atom_at == 'C':
            q_to_dist = -(group_charge - QM_res_charge)
        #q_to_dist = 0
        for x in zeroed:
            #print('q',df.loc[x,'q'])
            q_to_dist += df.loc[x, 'q']
        dist_df = dist_df[['q', 'X', 'Y', 'Z']]
        if round(q_to_dist, 5) != 0:
            q_to_dist = q_to_dist / len(dist_df.index.tolist())
            # distribute charge to atoms in dist_df dataframe
            dist_df['q'] = dist_df['q'] + q_to_dist
        output.write('Sum of residue group (should be an int):' + str(dist_df.sum()[0]) + '\n')
        MM_DZ_df = dist_df[['q', 'X', 'Y', 'Z']] 
     #   print(MM_DZ_df)
        for x in MM_DZ_df.values.tolist():
            for n,y in enumerate(x):
                if n < 3:
                    MM_for_array.append(str(y))
                else:
                    MM_for_array.append(str(y)+'\n')
    all_MM_atoms =  SEE_atoms(num_bonds_broken, with_HL) + with_HL['MM']
    all_MM_df = df.loc[df.index.isin(all_MM_atoms)]
    remain_MM_df = all_MM_df.drop(MM_used)
    ext_with_wat = remain_MM_df.index.tolist()
    for idx in remain_MM_df.index:
        if remain_MM_df.loc[idx, 'MOL2_RES'] == 'WAT' and remain_MM_df.loc[idx, 'MOL2_AT'] == 'OW':
            df_idx = remain_MM_df.loc[idx, 'MOL2_ID'] + 3.5
            ext_with_wat.append(df_idx)
    ext_wat_df = df.loc[df.index.isin(ext_with_wat)]
    remain_MM_df = ext_wat_df[['q', 'X', 'Y', 'Z']]
    for x in remain_MM_df.values.tolist():
        for n,y in enumerate(x):
            if n < 3:
                MM_for_array.append(str(y))
            else:
                MM_for_array.append(str(y)+'\n')
    #print(MM_used)
    mm_env = MM_for_array
    output.close()
    return mm_env


def Z2_atoms_zero(num_bonds_broken:int, with_HL:Dict[str,List[int]]) -> List[int]:
    """
    returns the M1  and M2 atoms that need to be zeroed

    Parameters
    ----------
    num_bonds_broken: int
        number of bonds broken that need to be capped
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region

    Returns
    -------
    MM_atoms: List[int]
        list of the MM atoms to zero
    """
    MM_atoms = []
    for bond in range(1,num_bonds_broken+1):
        if f'M1_{bond}' in with_HL.keys():
            MM_atoms += with_HL[f'M1_{bond}'] 
        if f'M2_{bond}' in with_HL.keys():
            MM_atoms += with_HL[f'M2_{bond}']
    return MM_atoms

def Z2_atoms_charge(num_bonds_broken:int, with_HL:Dict[str,List[int]]) -> List[str]:
    """
    returns the and M3 atoms 

    Parameters
    ----------
    num_bonds_broken: int
        number of bonds broken that need to be capped
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region

    Returns
    -------
    MM_atoms: List[int]
        list of the MM atoms to keep
    """
    MM_atoms = []
    for bond in range(1,num_bonds_broken+1):
        if f'M3_{bond}' in with_HL.keys():
            MM_atoms += with_HL[f'M3_{bond}']
    return MM_atoms

def Z2(df:pd.DataFrame, with_HL:Dict[str,List[int]], num_bonds_broken:int) -> List[str]:
    """
    zeroes the M1  and M2 atoms

    Parameters
    ----------
    df: pd.DataFrame
        pandas dataframe containing all atoms in the system
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region
    num_bonds_broken: int
        number of bonds broken that need to be capped

    Returns
    -------
    MM_for_array: List[str]
        list of the charges along with xyz coords
    """
    MM_for_array = []
    z2atoms0 = Z2_atoms_zero(num_bonds_broken, with_HL)
    z2atomsc = Z2_atoms_charge(num_bonds_broken, with_HL) + with_HL['MM']
    ext_with_wat = z2atoms0 + z2atomsc
    for idx in ext_with_wat:
        if df.loc[idx, 'MOL2_RES'] == 'WAT' and df.loc[idx, 'MOL2_AT'] == 'OW':
            df_idx = df.loc[idx, 'MOL2_ID'] + 3.5
            ext_with_wat.append(df_idx)
    ext_wat_df = df.loc[df.index.isin(ext_with_wat)]
    for idx in ext_wat_df.index:
        if idx in z2atoms0:
            ext_wat_df.at[idx, 'q'] = 0
    MM_df = ext_wat_df[['q', 'X', 'Y', 'Z']]
    for x in MM_df.values.tolist():
        for n,y in enumerate(x):
            if n < 3:
                MM_for_array.append(str(y))
            else:
                MM_for_array.append(str(y)+'\n')
    mm_env = MM_for_array
    return mm_env

def Z3_atoms_zero(num_bonds_broken:int, with_HL:Dict[str,List[int]]) -> List[int]:
    """
    returns the M1 and M2 and M3 atoms that need to be zeroed

    Parameters
    ----------
    num_bonds_broken: int
        number of bonds broken that need to be capped
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region

    Returns
    -------
    MM_atoms: List[int]
        list of the MM atoms to zero
    """
    MM_atoms = []
    for bond in range(1,num_bonds_broken+1):
        if f'M1_{bond}' in with_HL.keys():
            MM_atoms += with_HL[f'M1_{bond}'] 
        if f'M2_{bond}' in with_HL.keys():
            MM_atoms += with_HL[f'M2_{bond}'] 
        if f'M3_{bond}' in with_HL.keys():
            MM_atoms += with_HL[f'M3_{bond}']
    return MM_atoms

def Z3(with_HL:Dict[str,List[int]], num_bonds_broken:int) -> List[str]:
    """
    zeroes the M1 and M2  and M3 atoms

    Parameters
    ----------
    df: pd.DataFrame
        pandas dataframe containing all atoms in the system
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region
    num_bonds_broken: int
        number of bonds broken that need to be capped

    Returns
    -------
    MM_for_array: List[str]
        list of the charges along with xyz coords
    """
    MM_for_array = []
    z3atoms0 = Z3_atoms_zero(num_bonds_broken, with_HL)
    z3atomsc = with_HL['MM']
    ext_with_wat = z3atoms0 + z3atomsc
    for idx in ext_with_wat:
        if df.loc[idx, 'MOL2_RES'] == 'WAT' and df.loc[idx, 'MOL2_AT'] == 'OW':
            df_idx = df.loc[idx, 'MOL2_ID'] + 3.5
            ext_with_wat.append(df_idx)
    ext_wat_df = df.loc[df.index.isin(ext_with_wat)]
    for idx in ext_wat_df.index:
        if idx in z3atoms0:
            ext_wat_df.at[idx, 'q'] = 0
    MM_df = ext_wat_df[['q', 'X', 'Y', 'Z']]
    for x in MM_df.values.tolist():
        for n,y in enumerate(x):
            if n < 3:
                MM_for_array.append(str(y))
            else:
                MM_for_array.append(str(y)+'\n')
    mm_env = MM_for_array
    return mm_env

def get_charge_residue(MOL2_PATH:str, MM_line:str) -> Tuple[float, str]: # get charge for one MM PDB line
    """
    returns the charge of a specified MM atom based on coordinates

    Parameters
    ----------
    MOL2_PATH: str
        path to the mol2 file
    MM_line: str
        mm line to find in the mol2

    Returns
    -------
    charge, residue: Tuple[int,str]
        charge and residue of the MM line specified
    """
    with open(MOL2_PATH, 'r', encoding="iso-8859-1") as mol2file:
        lines = mol2file.readlines()
        for n,l in enumerate(lines):                                        # filter lines for atoms only
            if 'ATOM' in l:
                start = n
            elif 'HETATM' in l:
                start = n
            elif 'BOND' in l:
                end = n
        mollines = lines[start+1:end]
        for l in mollines:                                                  # matching lines in .mol2 file to MM XYZ coords bc
            mol2x = float("{:.3f}".format(float(l.split()[2])))             # indices between .mol2 and .pdb do not match
            mol2y = float("{:.3f}".format(float(l.split()[3])))
            mol2z = float("{:.3f}".format(float(l.split()[4])))
            if MM_line != None:
                for i in range(int(len(MM_line)/4)):                           # looping through each [atom,x,y,z]
                    mmx = float(MM_line[1+(4*i)])
                    mmy = float(MM_line[2+(4*i)])
                    mmz = float(MM_line[3+(4*i)])
                    if mmx == mol2x and mmy == mol2y and mmz == mol2z:
                        charge = l.split()[-2]
                        residue = l.split()[-3]+'_'+l.split()[-4]
        if 'charge' not in locals():
            charge = 0
            residue = ''
    return charge, residue

def redist_charges(num_bonds_broken:int, MM_for_array:List[str], MOL2_PATH:str, CAPPED_PDB_PATH:str, with_HL:Dict[str,List[int]]) -> List[str]:
    """
    RC: M1 charge is evenly redistributed to the midpoints of M1-M2 bonds

    Parameters
    ----------
    num_bonds_broken: int
        number of bonds broken that need to be capped
    MM_for_array: List[str]
        List of charges of each MM atom along with their XYZ coordinates
    MOL2_PATH: str
        path to mol2 file
    CAPPED_PDB_PATH: str
        path to capped complex pdb
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region

    Returns
    -------
    MM_for_array: List[str]
        list of the charges along with xyz coords
    """
    for bond in range(1,num_bonds_broken+1):
        # get coordinates of M1 and M2 atoms
        M1_atom = with_HL[f'M1_{bond}']
        M1atom_pdb = atoms_to_pdb_lines(CAPPED_PDB_PATH, M1_atom)
        M1atom_xyz = pdb_to_xyz(M1atom_pdb)
        M1atom_charge, M1atom_residue = get_charge_residue(MOL2_PATH, M1atom_xyz)  # this is for non-balanced scheme
        # M1atom_charge = q_total - q_i - q_QM
        M2_atoms = with_HL[f'M2_{bond}']
        num_M2 = len(M2_atoms)
        redist_charge = float(M1atom_charge) / float(num_M2)
        for atom in M2_atoms:
            midpoint_xyz = []
            M2atom_pdb = atoms_to_pdb_lines(CAPPED_PDB_PATH, [atom])
            M2atom_xyz = pdb_to_xyz(M2atom_pdb)
            M2atom_charge, M2atom_residue = get_charge_residue(MOL2_PATH, M2atom_xyz)
            # appending M2 atoms to array with new charge for RCD only
            if charge_method == 'RCD':
                rcd_m2_charge = float(M2atom_charge) - float(redist_charge)
                MM_for_array.append(f'{float(rcd_m2_charge):.3f}')
                MM_for_array.append(f'{float(M2atom_xyz[1]):.3f}')
                MM_for_array.append(f'{float(M2atom_xyz[2]):.3f}')
                MM_for_array.append(f'{float(M2atom_xyz[3]):.3f}\n')
            if charge_method == 'RC':
                MM_for_array.append(str(redist_charge))
            elif charge_method == 'RCD':
                MM_for_array.append(str(2*redist_charge)) # double redistributed charge for RCD scheme
            # calculate midpoints of M1 and M2
            for n,x in enumerate(M2atom_xyz[1:]):
                coord = (float(M1atom_xyz[n+1]) + float(x))/2
                midpoint_xyz.append(f'{coord:.3f}')
                if n == 0 or n == 1:
                    MM_for_array.append(f'{coord:.3f}')
                else:
                    MM_for_array.append(f'{coord:.3f}\n')
    return MM_for_array


def bal_redist_charges(num_bonds_broken:int, MM_for_array:List[str], MM_atoms:List[str], charge_method:str, df:pd.DataFrame, with_HL:Dict[str,List[int]]) -> List[str]:
    """
    TODO Insert description

    Parameters
    ----------
    num_bonds_broken: int
        number of bonds broken that need to be capped
    MM_for_array: List[str]
        List of charges of each MM atom along with their XYZ coordinates
    MM_atoms: MM_atoms: List[str]
       List of MM atoms 
    charge_method: str
        redistribution scheme
    df: pd.DataFrame
        pandas dataframe containing all atoms in the system
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region

    Returns
    -------
    MM_for_array: List[str]
        list of the charges along with xyz coords
    """
    for bond in range(1,num_bonds_broken+1):
        # get coordinates of M1 and M2 atoms
        M1_atom = with_HL[f'M1_{bond}']
        M1atom_xyz = [df.at[M1_atom[0], 'X'], df.at[M1_atom[0], 'Y'], df.at[M1_atom[0],'Z']]
        M1atom_charge = df.at[M1_atom[0], 'q']
        M1atom_residue = df.at[M1_atom[0], 'PDB_RES']
        if df.at[M1_atom[0], 'PDB_AT'] == 'CA':     
            # get original integer charge of residue
            MM_group_charges = np.where(df['PDB_RES'] == str(M1atom_residue), df[['q']].sum(axis=1),0)
            MM_group_charge = np.sum(MM_group_charges)
            #print('MM group charge (int):', MM_group_charge)
            # sum the charges on atoms in that residue and also in MM
            MMdf = df.loc[df.index.isin(MM_atoms)]
            MM_resi_charges = np.where(MMdf['PDB_RES'] == str(M1atom_residue), MMdf[['q']].sum(axis=1),0)
            MM_resi_charge = np.sum(MM_resi_charges)
            #print('MM_resi_charge:',MM_resi_charge)
            q_purple = MM_group_charge - MM_resi_charge
            #print('q_purple:', q_purple)
            M1_bal_charge = M1atom_charge + q_purple
            #print('New charge:', M1_bal_charge + MM_resi_charge)
        elif df.at[M1_atom[0], 'PDB_AT'] == 'C':
            #print('M1 atom type is C')
            # get original integer charge of residue
            QM_group_charges = np.where(df['PDB_RES'] == str(M1atom_residue), df[['q']].sum(axis=1),0)
            QM_group_charge = np.sum(QM_group_charges)
            #print('QM group charge (int):', QM_group_charge)
            # sum the charges on atoms in that residue and also in MM
            MMdf = df.loc[df.index.isin(MM_atoms)]
            MM_resi_charges = np.where(MMdf['PDB_RES'] == str(M1atom_residue), MMdf[['q']].sum(axis=1),0)
            q_purple = np.sum(MM_resi_charges)
            #print('q_purple:',q_purple)
            M1_bal_charge = M1atom_charge - q_purple
        #print('M1_bal_charge:', M1_bal_charge)
        M2_atoms = with_HL[f'M2_{bond}']
        num_M2 = len(M2_atoms)
        redist_charge = float(M1_bal_charge) / float(num_M2)
        #print('Redist charge:', redist_charge)
        for atom in M2_atoms:
            midpoint_xyz = []
            M2atom_charge = df.at[atom, 'q']
            M2atom_xyz = [df.at[atom, 'X'], df.at[atom, 'Y'], df.at[atom,'Z']]
            # appending M2 atoms to array with new charge for RCD only
            if charge_method == 'BRCD':
                rcd_m2_charge = float(M2atom_charge) - float(redist_charge)
                MM_for_array.append(f'{float(rcd_m2_charge):.6f}')
                MM_for_array.append(f'{float(M2atom_xyz[0]):.3f}')
                MM_for_array.append(f'{float(M2atom_xyz[1]):.3f}')
                MM_for_array.append(f'{float(M2atom_xyz[2]):.3f}\n')
            if charge_method == 'BRC2':
                rcd_m2_charge = float(M2atom_charge) + float(redist_charge)
                MM_for_array.append(f'{float(rcd_m2_charge):.6f}')
                MM_for_array.append(f'{float(M2atom_xyz[0]):.3f}')
                MM_for_array.append(f'{float(M2atom_xyz[1]):.3f}')
                MM_for_array.append(f'{float(M2atom_xyz[2]):.3f}\n')
            if charge_method == 'BRC':
                MM_for_array.append(f'{float(redist_charge):.6f}')
            elif charge_method == 'BRCD':
                MM_for_array.append(f'{float(2*redist_charge):.6f}') # double redistributed charge for RCD scheme
            # calculate midpoints of M1 and M2
            if charge_method == 'BRC' or charge_method == 'BRCD':
                for n,x in enumerate(M2atom_xyz):
                    coord = (float(M1atom_xyz[n]) + float(x))/2
                    midpoint_xyz.append(f'{coord:.3f}')
                    if n == 0 or n == 1:
                        MM_for_array.append(f'{coord:.3f}')
                    else:
                        MM_for_array.append(f'{coord:.3f}\n')
    return MM_for_array

def balanced_RC(charge_method:str, df:pd.DataFrame, with_HL:Dict[str,List[int]], num_bonds_broken:int) -> List[str]:
    """
    TODO Insert description

    Parameters
    ----------
    charge_method: str
        redistribution scheme
    df: pd.DataFrame
        pandas dataframe containing all atoms in the system
    with_HL: Dict[str, List[int]]
        Dictionary with key of the region and value with a list of atoms in that region
    num_bonds_broken: int
        number of bonds broken that need to be capped

    Returns
    -------
    MM_for_array: List[str]
        list of the charges along with xyz coords
    """
    #print(df)
    if charge_method == 'BRC':
        BRC_ext = Z1_atoms_charge(num_bonds_broken, with_HL) + with_HL['MM']
    elif charge_method == 'BRCD' or charge_method == 'BRC2':
        BRC_ext = Z2_atoms_charge(num_bonds_broken, with_HL) + with_HL['MM']
    elif charge_method == 'BRC3':
        BRC_ext = Z3_atoms_charge(num_bonds_broken) + with_HL['MM']
    ext_df = df.loc[df.index.isin(BRC_ext)]
    # add fourth water point
    ext_with_wat = BRC_ext
    for idx in ext_df.index:
        if ext_df.loc[idx, 'MOL2_RES'] == 'WAT' and ext_df.loc[idx, 'MOL2_AT'] == 'OW':
            df_idx = ext_df.loc[idx, 'MOL2_ID'] + 3.5
            ext_with_wat.append(df_idx)
    ext_wat_df = df.loc[df.index.isin(ext_with_wat)]
    ext_wat_df.to_csv('MMdf.csv')
    ext_wat_df = ext_wat_df[['q', 'X', 'Y', 'Z']]
    MM_for_array = []
    for x in ext_wat_df.values.tolist():
        for n,y in enumerate(x):
            if n < 3:
                MM_for_array.append(str(y))
            else:
                MM_for_array.append(str(y)+'\n')
    MM_atoms =  SEE_atoms(num_bonds_broken, with_HL) + with_HL['MM']
    mm_env = bal_redist_charges(num_bonds_broken, MM_for_array, MM_atoms, charge_method, df, with_HL)
    #print(mm_env)
    return mm_env

#TODO do we want to keep this method
def write_extern_xyz(filepath:str, mm_env:List[str]) -> None:
    """
    Given a list of point charges and the coordinates to place those charges, writes an xyz 
    
    Parameters
    ----------
    filepath: str
        filepath to write the external xyz
    MM_env: List[str]
        list of all the MM atoms with the point charges and x, y, z coordinates

    Returns
    -------
    None
    """
    with open(filepath, 'w') as wfile:
        wfile.write(' '.join(mm_env)[:-1])

#def write_psi4_file(PSI4_FILE_PATH):


def write_QM(charge_method:str, c_ligand:str, basis_set:str, method:str, PSI4_FILE_PATH:str, do_fsapt: bool = None) -> None:
    """
    creates psi4 python input files for different charge schemes:
    SEE: no changes to point charges (M1 charge does not change)
    Z1: M1 charge is zeroed-out
    Z2: M1 and M2 charges are zeroed-out
    Z3: M1, M2, and M3 charges are zeroed-out
    RC: M1 charge is evenly redistributed to the midpoints of M1-M2 bonds
    RCD: same as RC, but M1-M2 bond dipoles are preserved

    Parameters
    ----------
    charge_method: str
        redistribution scheme
    c_ligand: str
        charge of the ligand
    basis_set: str
        basis set for the QM computation
    method: str
        method for QM energy 
    PSI4_FILE_PATH: str
        name for created psi4 file
    do_fsapt: boolean or None
        if fsapt needs to be turned off in the psi4 input file this option will be False


    Returns
    -------
    None
    """
    CAPPED_PDB_PATH = 'CAPPED-prot_autocap_fixed.pdb'
    LIGAND_PDB_PATH = 'ligand.pdb'
    CAPPED_QM_PDB_PATH = 'CAPPED_qm.pdb'
    MOL2_PATH = 'prot_autocap_fixed.mol2'
    df = pd.read_csv('dataframe.csv', index_col=['CX_PDB_ID'])
    with open('with_HL.dat', 'r') as dictfile:
        with_HL = json.load(dictfile)
    num_bonds_broken = int(list(with_HL.keys())[-1].split('_')[-1])
    c_QM = 0
    for x in with_HL['QM']:
        if '+' in df.at[x,'AT_LABEL']:
            c_QM += int(df.at[x, 'AT_LABEL'][-2])
        elif '-' in df.at[x, 'AT_LABEL']:
            c_QM -= int(df.at[x, 'AT_LABEL'][-2])
    c_QM = str(int(c_QM))

    #print('num_bonds_broken:', num_bonds_broken)
    inp_filename = PSI4_FILE_PATH.split('/')[-1]
    MM_for_array = []
    # get QM XYZ
    QM_atoms = with_HL['QM']
    for bond in range(1,num_bonds_broken+1):
        QM_atoms += with_HL[f'Q1_{bond}']
        QM_atoms += with_HL[f'HL_{bond}']
    QM_atoms.sort()
    qm_pdb_pro_lines = atoms_to_pdb_lines(CAPPED_PDB_PATH, QM_atoms)
    write_pdb(CAPPED_PDB_PATH, QM_atoms, CAPPED_QM_PDB_PATH, LIGAND_PDB_PATH)
    qm_pro = pdb_to_xyz(qm_pdb_pro_lines)
    qm_pdb_lig_lines = qm_ligand(LIGAND_PDB_PATH)
    qm_lig = pdb_to_xyz(qm_pdb_lig_lines)
    # get MM array of frontier region
    if charge_method == 'SEE':
        mm_env = SEE(MOL2_PATH, CAPPED_PDB_PATH, with_HL, num_bonds_broken) #appends appropriate M1, M2, M3 XYZ lines to MM_for_array
    elif charge_method == 'Z1':
          mm_env = Z1(df, with_HL, num_bonds_broken)
    elif charge_method == 'DZ1':
          mm_env = DZn(1, df, with_HL, num_bonds_broken)
    elif charge_method == 'Z2':
          mm_env = Z2(df, with_HL, num_bonds_broken)
    elif charge_method == 'DZ2':
          mm_env = DZn(2, df, with_HL, num_bonds_broken)
    elif charge_method == 'Z3':
          mm_env = Z3(with_HL, num_bonds_broken)
    elif charge_method == 'DZ3':
          mm_env = DZn(3, df, with_HL, num_bonds_broken)
    elif charge_method == 'BRC':
          mm_env = balanced_RC(charge_method, df, with_HL, num_bonds_broken)
    elif charge_method == 'BRC2':
          mm_env = balanced_RC(charge_method, df, with_HL, num_bonds_broken)
    elif charge_method == 'BRCD':
          mm_env = balanced_RC(charge_method, df, with_HL, num_bonds_broken)
    else:
        #we check for this in the input file so this else shouldn't ever be reached
        print('incorrect charge scheme')
        sys.exit()
    
    #write_extern_xyz(f'{charge_method}/extern.xyz', mm_env)
    #dump_pkl()
    #write_psi4_file(PSI4_FILE_PATH, lig_charge)
    with open(PSI4_FILE_PATH, 'a') as inpfile:
        inpfile.write("""
import psi4
import numpy as np
import qcelemental as qcel 
import time

start = time.time()

psi4.set_memory('70 GB')
psi4.core.set_num_threads(10)

psi4.core.set_output_file('""" + f"{inp_filename[:-3]}.out'" + """, False)\n


dimer =psi4.geometry('''\n""" + c_ligand + ' 1\n'
+ ' '.join(qm_lig) + '--\n' + c_QM + ' 1\n '
+ ' '.join(qm_pro) + 
"""units angstrom
symmetry c1
no_com
no_reorient
''')\n
Chargefield_B = np.array([\n"""+ 
','.join(mm_env)[:-1] +
"""]).reshape((-1,4))\n""" +
"""Chargefield_B[:,[1,2,3]] /= qcel.constants.bohr2angstroms\n
psi4.set_options({
'basis': '""" + basis_set +"""',\n"""+
"""'freeze_core': 'True',
'scf_type': 'df',
'mp2_type': 'df'\n""")
        if do_fsapt == False:
            inpfile.write("'do_fsapt': 'False'\n")
        inpfile.write("""})\n
e = psi4.energy('"""+method+"""', external_potentials={'B':Chargefield_B})\n

end=time.time()
wall_time = '{:.2f}'.format(float(end-start))
with open ('"""+inp_filename[:-3]+""".out', 'a') as output:
    output.write(f'Wall time: {wall_time} seconds')""")
'''
def dump_pkl():
    results = {'external':[], 'ligand': [], 'protein': []}
    results['external'] = mm_env 
    results['ligand'] = c_ligand + ' 1\n'+' '.join(qm_lig) 
    results['protein'] = c_QM + ' 1\n' +' '.join(qm_pro)
    with open('results.dat','w+') as out:
        out.write(json.dumps(results))
'''

def write_input(inputfile, psi4file):
    """
    Writes the Sparcle-QC input file to the top of the psi4file

    Parameters
    ----------
    inputfile: str
        name of the Sparcle-QC input file
    psi4file: str
        name of created psi4 input file

    Returns
    -------
    None
    """
    with open(inputfile) as inp:
        with open(psi4file, 'w') as psi4file:
            psi4file.write('"""\nThis Psi4 file was created using Sparcle-QC with the following specifications:\n')
            for line in inp:
                psi4file.write(line)
            psi4file.write('"""\n\n')


def check_QM_file(psi4file: str) -> None:

    """
    Checks that the charge of the QM region in the created input file is integer
    prints the charge, along with the number of atoms in the QM and MM region

    Parameters
    ----------
    psi4file: str
        name of created psi4 input file

    Returns
    -------
    None
    """
    out = open(glob('*.out')[0], 'a')

    with open(psi4file, 'r') as pfile:
        lines = pfile.readlines()
        extern_idx = []
        dimer_idx = []
        prot_idx = []
        for n,l in enumerate(lines):
            if 'Chargefield' in l:
                extern_idx.append(n)
            elif 'dimer' in l:
                dimer_idx.append(n)
            elif '--' in l:
                prot_idx.append(n)
            elif 'unit' in l:
                dimer_idx.append(n)
        array = lines[extern_idx[0]+1:extern_idx[1]]
        charge = float(array[0].split(',')[0])
        num_atoms = dimer_idx[1] - dimer_idx[0] - 4 + 1
        num_prot_atoms = dimer_idx[1] - prot_idx[0] - 2
        for l in array[1:]:
            charge += float(l.split(',')[1])
            num_atoms += 1
        
        out.write('\n----------------------------------------------------------------------------------------------------\n')
        out.write('check_psi4_file'.center(100)+'\n')
        out.write('----------------------------------------------------------------------------------------------------\n')

        out.write(f'Psi4 filename: {psi4file}\n')
        out.write(f'Total charge of MM region: {charge:.2f}\n')
        out.write(f'Number of point charges in MM region: {num_atoms}\n')
        out.write(f'Number of non-ligand atoms in QM region: {num_prot_atoms}')
    
    out.close()
