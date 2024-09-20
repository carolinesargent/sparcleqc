import json
import pandas as pd
from typing import Dict, List

def get_M2_N_resn(bond:str, df:pd.DataFrame, no_HL:Dict[str, List[int]]) -> str:
    """
    Given a bond number, this function returns the resname that bond is apart of, indicated by the N atom

    Parameters
    ----------
    bond: str
        number of bond in HL to find resname for
    df: pd.DataFrame
        pandas dataframe to search for the resname in
    no_HL: Dict[str,List[int]]
        Dictionary indicating the atoms in each region

    Returns
    -------
    M2_N_resn: str
        resname for the corresponding M2 atom
    """
    M2atoms = no_HL[f'M2_{bond}']
    for atom in M2atoms:
        if df.loc[atom, 'PDB_AT'] == 'N':
            M2_N_resn = df.loc[atom, 'PDB_RES']
    return M2_N_resn

def find_M3(bond:str, M2resn:str, df:pd.DataFrame, no_HL:Dict[str,List[int]]) -> int:
    """
    Given a bond number, this function returns an atom number if there is an M3 atom for that bond that is a different resname than the M2resn

    Parameters
    ----------
    bond: str
        number bond to find M3 for
    M2resn: str
        name of residue for the M2 atom
    df: pd.DataFrame
        pandas dataframe to search for the resname in
    no_HL: Dict[str,List[int]]
        Dictionary indicated the atoms in each region

    Returns
    -------
    M3_to_move: int or None
        atom number for the M3 atom if it has a different residue name than the M2 atom
    """
    M3atoms = no_HL[f'M3_{bond}']
    for atom in M3atoms:
        if df.loc[atom, 'PDB_RES'] != M2resn:
            M3_to_move = atom
            return M3_to_move

def get_key(val:str, no_HL:Dict[str,List[int]]) -> str:
    """
    Given a value in the dictionary defining the atoms in each region, this function returns the key that corresponds to the given unique value

    Parameters
    ----------
    val: str
        value to search for in the dictionary
    no_HL: Dict[str,List[int]]
        Dictionary indicated the atoms in each region

    Returns
    -------
    key: str
        key in with_HL that the value is a part of
    """
    for key, values in no_HL.items():
        for value in values:
            if int(val) == value:
                return key

def move_ids(atom_ids: List[int], no_HL:Dict[str,List[int]]) -> Dict[str,List[int]]:
    """
    Moves the specified M3 atoms into the MM region of the dictionary

    Parameters
    ----------
    atom_ids: List[int]
        atoms to move from M3 to M
    no_HL: Dict[str,List[int]]
        Dictionary indicating the atoms in each region

    Returns
    -------
    no_HL: Dict[str,List[int]]
        updated dictionary with the M3 atoms moved to the M region
    """
    for atom_id in atom_ids:
        old_k = get_key(atom_id, no_HL)
        if 'M3' in old_k:
            no_HL['MM'].append(atom_id)
            no_HL[old_k].remove(atom_id)
    return no_HL


def move_m3s() -> None:    
    """
    This function moves atoms from the M3 region to the MM region in dictionary.dat
    only if the M3 atom has a different residue than the rest of the MM residue. 
    This residue is identified by the M3 N atom. 

    Parameters
    ----------
    None

    Returns
    -------
    None
    """
    with open('pre-dictionary.dat', 'r') as dictfile:
        no_HL = json.load(dictfile)
    num_broken_bonds = int(list(no_HL.keys())[-1].split('_')[-1])
    df = pd.read_csv('dataframe.csv', index_col=['CX_PDB_ID'])
    to_move = []
    for bond in range(1,num_broken_bonds+1):
        M2resn = get_M2_N_resn(bond, df, no_HL)
        atom_to_move = find_M3(bond, M2resn, df, no_HL)
        if atom_to_move is not None:
            to_move.append(find_M3(bond, M2resn, df, no_HL))
    new_dict = move_ids(to_move, no_HL)
    with open('dictionary.dat', 'w+') as dictfile:
        dictfile.write(json.dumps(new_dict))
