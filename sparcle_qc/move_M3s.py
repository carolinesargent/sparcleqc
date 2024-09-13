'''
This script moves atoms from the M3 region to the MM region in dictionary.dat
only if the M3 atom has a different residue than the rest of the MM residue. 
This residue is identified by the M3 N atom. 
'''
import json
import pandas as pd


def get_M2_N_resn(bond, df, no_HL):
    M2atoms = no_HL[f'M2_{bond}']
    for atom in M2atoms:
        if df.loc[atom, 'PDB_AT'] == 'N':
            M2_N_resn = df.loc[atom, 'PDB_RES']
    return M2_N_resn

def find_M3(bond, M2resn, df, no_HL):
    M3atoms = no_HL[f'M3_{bond}']
    for atom in M3atoms:
        if df.loc[atom, 'PDB_RES'] != M2resn:
            M3_to_move = atom
            return M3_to_move

def get_key(val, no_HL):
    for key, values in no_HL.items():
        for value in values:
            if int(val) == value:
                return key

def move_ids(atom_ids, no_HL):
    for atom_id in atom_ids:
        old_k = get_key(atom_id, no_HL)
        if 'M3' in old_k:
            no_HL['MM'].append(atom_id)
            no_HL[old_k].remove(atom_id)
    return no_HL


def move_m3s():    
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
