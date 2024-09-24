import pandas as pd
from typing import Dict

def read_pdb(PDB_PATH: str, d:Dict) -> pd.DataFrame:
    """ 
    Populates a pandas dataframe with the information in the pdb
    (x_coord, y_coord, z_coordi, resname, atom type) indexed by atom id

    Parameters
    ----------
    PDB_PATH: str
        path to protein pdb to populate the dataframe
    d: Dict
        Dictionary that has one key, PDB_ID, with the value as an empty list

    Returns
    -------
    df: pd.DataFrame
        dataframe containing the information from the PDB

    """
    with open(PDB_PATH, 'r') as pdbfile:
        lines = pdbfile.readlines()
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
                d['PDB_ID'].append(l[6:11].strip())
        df = pd.DataFrame(d)
        df = df.set_index('PDB_ID')
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
                df.loc[l[6:11].strip(), 'PDB_AT'] = l[11:16].strip()

                x_coord = l[29:38].strip()
                y_coord = l[38:46].strip()
                z_coord = l[46:54].strip()

                df.loc[l[6:11].strip(), 'PDB_RES'] = l[16:20].strip()+'_'+l[22:26].strip()
                df.loc[l[6:11].strip(), 'X'] = float(x_coord)
                df.loc[l[6:11].strip(), 'Y'] = float(y_coord)
                df.loc[l[6:11].strip(), 'Z'] = float(z_coord)
                df.loc[l[6:11].strip(), 'AT_LABEL'] = l[66:87].strip() 
    return df

def read_mol2(MOL2_PATH:str, m:Dict) -> pd.DataFrame:
    """ 
    Populates a pandas dataframe with the information in the mol2
    (x_coord, y_coord, z_coordi, resname, atom type) indexed by atom id

    Parameters
    ----------
    MOL2_PATH: str
        path to protein mol2 to populate the dataframe
    d: Dict
        Dictionary that has one key, MOL2_ID, with the value as an empty list

    Returns
    -------
    df: pd.DataFrame
        dataframe containing the information from the MOL2

    """
    with open(MOL2_PATH, 'r') as molfile:
        lines = molfile.readlines()
        for l in lines:
            if len(l.split()) > 7:
                m['MOL2_ID'].append(l.split()[0])
        df = pd.DataFrame(m)
        df = df.set_index('MOL2_ID')
        for l in lines:
            if len(l.split()) > 7:
                df.loc[l.split()[0], 'MOL2_AT'] = l.split()[5]
                df.loc[l.split()[0], 'MOL2_RES'] = l.split()[-3]
                x_mol = float("{:.3f}".format(float(l.split()[2])))
                y_mol = float("{:.3f}".format(float(l.split()[3])))
                z_mol = float("{:.3f}".format(float(l.split()[4])))
                df.loc[l.split()[0], 'X'] = x_mol
                df.loc[l.split()[0], 'Y'] = y_mol
                df.loc[l.split()[0], 'Z'] = z_mol
                df.loc[l.split()[0], 'q'] = l.split()[-2]
    return df        

def combine(pdb_df:pd.DataFrame, mol2_df:pd.DataFrame) -> pd.DataFrame:
    """ 
    combines information from a pdb dataframe with a
    mol2 dataframe based on matching the coordinates

    Parameters
    ----------
    pdb_df: pd.DataFrame
        dataframe containing information from the protein pdb
    mol2_df: pd.DataFrame
        dataframe containing information from the protein mol2

    Returns
    -------
    df: pd.DataFrame
        dataframe containing information combined from the pdb_df and the mol2_df with one entry per atom

    """
    for idx in pdb_df.index:
        mol2_idx_x = mol2_df.loc[mol2_df['X'] == pdb_df.loc[idx, 'X']].index.tolist()
        mol2_idx_y = mol2_df.loc[mol2_df['Y'] == pdb_df.loc[idx, 'Y']].index.tolist()
        mol2_idx_z = mol2_df.loc[mol2_df['Z'] == pdb_df.loc[idx, 'Z']].index.tolist()
        mol2_idx_xyz = mol2_idx_x + mol2_idx_y + mol2_idx_z
        for x in mol2_idx_xyz:
            if mol2_idx_xyz.count(x) == 3:
                mol2_idx = x
                pdb_df.loc[idx,'MOL2_ID'] = mol2_idx
                pdb_df.loc[idx,'MOL2_AT'] = mol2_df.loc[mol2_idx,'MOL2_AT']
                pdb_df.loc[idx,'MOL2_RES'] = mol2_df.loc[mol2_idx,'MOL2_RES']
                pdb_df.loc[idx,'q'] = mol2_df.loc[mol2_idx,'q']
        # insert fourth point of water, if present in mol2
        if 'HOH' in pdb_df.loc[idx, 'PDB_RES'] and pdb_df.loc[idx, 'PDB_AT'] == 'O':
            EPW_idx = str(int(mol2_idx)+3)
            if EPW_idx in mol2_df.index and 'EP' in mol2_df.loc[EPW_idx, 'MOL2_AT']:
                pdb_df.loc[str(float(EPW_idx)+.5), 'MOL2_AT'] = mol2_df.loc[EPW_idx, 'MOL2_AT']
                pdb_df.loc[str(float(EPW_idx)+.5), 'MOL2_RES'] = mol2_df.loc[EPW_idx, 'MOL2_RES']
                pdb_df.loc[str(float(EPW_idx)+.5), 'MOL2_ID'] = EPW_idx
                pdb_df.loc[str(float(EPW_idx)+.5), 'q'] = mol2_df.loc[EPW_idx, 'q']
                pdb_df.loc[str(float(EPW_idx)+.5), 'X'] = mol2_df.loc[EPW_idx, 'X']
                pdb_df.loc[str(float(EPW_idx)+.5), 'Y'] = mol2_df.loc[EPW_idx, 'Y']
                pdb_df.loc[str(float(EPW_idx)+.5), 'Z'] = mol2_df.loc[EPW_idx, 'Z']
                pdb_df.loc[str(float(EPW_idx)+.5), 'PDB_RES'] = pdb_df.loc[idx, 'PDB_RES']
    return pdb_df

def read_cx_pdb(CX_PDB_PATH:str, d:Dict) -> pd.DataFrame:
    """ 
    Populates a pandas dataframe with the information in the pdb
    (x_coord, y_coord, z_coordi, resname, atom type) indexed by atom id

    Parameters
    ----------
    CX_PDB_PATH: str
        path to complex pdb to populate the dataframe
    d: Dict
        Dictionary that has one key, CX_PDB_ID, with the value as an empty list

    Returns
    -------
    df: pd.DataFrame
        dataframe containing the information from the PDB

    """
    with open(CX_PDB_PATH, 'r') as pdbfile:
        lines = pdbfile.readlines()
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
                d['CX_PDB_ID'].append(l[6:11].strip())
        df = pd.DataFrame(d)
        df = df.set_index('CX_PDB_ID')
        for l in lines:
            if l[0:6].strip() == 'ATOM' or l[0:6].strip() == 'HETATM':
                x_coord = l[29:38].strip()
                y_coord = l[38:46].strip()
                z_coord = l[46:54].strip()
                df.loc[l[6:11].strip(), 'X'] = float(x_coord)
                df.loc[l[6:11].strip(), 'Y'] = float(y_coord)
                df.loc[l[6:11].strip(), 'Z'] = float(z_coord)
    return df

def combine2(pdb_df:pd.DataFrame, cx_pdb_df:pd.DataFrame) -> pd.DataFrame:
    """
    combines information from a dataframe from the protein pdb with a
    dataframe from the cx pdb based on matching the coordinates

    Parameters
    ----------
    pdb_df: pd.DataFrame
        dataframe containing information from the protein pdb
    cx_pdb_df: pd.DataFrame
        dataframe containing information from the cx pdb

    Returns
    -------
    df: pd.DataFrame
        dataframe containing information combined from the pdb_df and the cx_df with one entry per atom

    """
    for idx in pdb_df.index:
        cx_idx_x = cx_pdb_df.loc[cx_pdb_df['X'] == pdb_df.loc[idx, 'X']].index.tolist()
        cx_idx_y = cx_pdb_df.loc[cx_pdb_df['Y'] == pdb_df.loc[idx, 'Y']].index.tolist()
        cx_idx_z = cx_pdb_df.loc[cx_pdb_df['Z'] == pdb_df.loc[idx, 'Z']].index.tolist()
        cx_idx_xyz = cx_idx_x + cx_idx_y + cx_idx_z
        for x in cx_idx_xyz:
            if cx_idx_xyz.count(x) == 3:
                cx_idx = x
                pdb_df.loc[idx,'CX_PDB_ID'] = cx_idx
        if pdb_df.loc[idx, 'MOL2_AT'] == 'EP':
            pdb_df.loc[idx,'CX_PDB_ID'] = idx
    return pdb_df

def change_water_charges(df: pd.DataFrame, o: str, h: str, ep:str = None) -> pd.DataFrame:
    """
    Changes charges of the water atoms (and possibly EP) depending on the user specified inputs

    Parameters
    ----------
    df: pd.DataFrame
        dataframe containing information from the mol2s and pdbs
    o: str
        desired oxygen charge
    h: str
        desired hydrogen charge
    ep: str
        desired extra point charge

    Returns
    -------
    df: pd.DataFrame
        original dataframe with water charges updated accordingly 

    """
    for idx in df.index:
        resi = df.loc[idx, 'PDB_RES']
        if 'HOH' in resi or 'WAT' in resi:
            atom = df.loc[idx, 'MOL2_AT']
            if 'O' in atom:
                df.loc[idx, 'q'] = float(o)
            elif 'H' in atom:
                df.loc[idx, 'q'] = float(h)
            elif 'EP' in atom:
                if ep != None:
                    df.loc[idx, 'q'] = float(ep)
    return df

def combine_data(o_charge:str = None, h_charge:str = None, ep_charge:str = None) -> None:
    """
    Creates a dataframe from the provided protein pdb, complex pdb,
    and mol2 information and then updates it with the desired water
    atom charges

    Parameters
    ----------
    o: str
        desired oxygen charge
    h: str
        desired hydrogen charge
    ep: str
        desired extra point charge

    Returns
    -------

    """
    PDB_PATH = 'prot_autocap_fixed.pdb'
    MOL2_PATH = 'prot_autocap_fixed.mol2'
    CX_PDB_PATH = 'cx_autocap_fixed.pdb'
    p = {'PDB_ID':[]}
    m = {'MOL2_ID':[]}
    c = {'CX_PDB_ID':[]}
    pdb_info = read_pdb(PDB_PATH, p)
    mol2_info = read_mol2(MOL2_PATH, m)
    cx_pdb_info = read_cx_pdb(CX_PDB_PATH, c)
    combined = combine(pdb_info, mol2_info)
    combined2 = combine2(combined, cx_pdb_info)
    if o_charge!= None and h_charge != None and ep_charge !=None:
        final = change_water_charges(combined2, o_charge, h_charge, ep_charge)
    elif o_charge != None and h_charge != None:
        final = change_water_charges(combined2, o_charge , h_charge)
    else:
        final = combined2
    final.to_csv('dataframe.csv')
