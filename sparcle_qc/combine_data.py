import pandas as pd
import json
import os
import sys


def read_pdb(PDB_PATH, d):
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

def read_mol2(MOL2_PATH, m): 
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

def combine(pdb_df, mol2_df):
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
        # insert fourth point of water
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

def read_cx_pdb(CX_PDB_PATH, d):
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

def combine2(pdb_df, cx_pdb_df):
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

def check_charges(combined_df, mol2_info):
    # Get indicies of Nan or zero charges
    df = combined_df[['q']]
    d = df.apply(lambda s: pd.to_numeric(s, errors="coerce"))
    m = d.eq(0) | d.isna()
    s = m.stack()
    indices = s[s].index.tolist()
    #print(indices)
    # For Hs of ACE or NME, correct their charges. Else, print residue.
    for idx, q in indices:
        pdb_at = combined_df.loc[idx, 'PDB_AT']
        mol2_at = combined_df.loc[idx, 'MOL2_AT']
        pdb_res = combined_df.loc[idx, 'PDB_RES'][:3]
        
        '''
        if 'ACE' in pdb_res and pdb_at == '1HH3':
            combined_df.loc[idx, 'q'] = 0.1123
        if 'ACE' in pdb_res and pdb_at == '2HH3':
            combined_df.loc[idx, 'q'] = 0.1123
        if 'ACE' in pdb_res and pdb_at == '3HH3':
            combined_df.loc[idx, 'q'] = 0.1123
        if 'NME' in pdb_res and pdb_at == '1HH3':
            combined_df.loc[idx, 'q'] = 0.0976
        if 'NME' in pdb_res and pdb_at == '2HH3':
            combined_df.loc[idx, 'q'] = 0.0976
        if 'NME' in pdb_res and pdb_at == '3HH3':
            combined_df.loc[idx, 'q'] = 0.0976
        if 'NME' in pdb_res and pdb_at == 'CH3':
            combined_df.loc[idx, 'q'] = -0.1490
        '''
        #if pdb_res != 'NME' and pdb_res != 'ACE' and pdb_res != 'HOH':
            #print('For the residue(s) below, check charges between the mol2 file and dataframe.csv. Dataframe.csv likely needs updating.')
            #print(combined_df.loc[idx, 'PDB_RES'])
    return combined_df

def change_water_charges(df, o, h, ep=None):
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

def combine_data(o_charge=None, h_charge = None, ep_charge = None):
    PDB_PATH = 'prot_autocap_fixed.pdb'
    MOL2_PATH = 'prot_autocap_fixed.mol2'
    CX_PDB_PATH = 'cx_autocap_fixed.pdb'
    p = {'PDB_ID':[]}
    m = {'MOL2_ID':[]}
    c = {'CX_PDB_ID':[]}
    pdb_info = read_pdb(PDB_PATH, p)
    mol2_info = read_mol2(MOL2_PATH, m)
    cx_pdb_info = read_cx_pdb(CX_PDB_PATH, c)
    #print(pdb_info.head(20))
    #print(mol2_info.head(20))
    #print(cx_pdb_info.head(20))
    combined = combine(pdb_info, mol2_info)
    #print(combined.head(20))
    combined2 = combine2(combined, cx_pdb_info)
    if o_charge!= None and h_charge != None and ep_charge !=None:
        final = change_water_charges(combined2, o_charge, h_charge, ep_charge)
    elif o_charge != None and h_charge != None:
        final = change_water_charges(combined2, o_charge , h_charge)
    else:
        final = combined2
    final.to_csv('dataframe.csv')
