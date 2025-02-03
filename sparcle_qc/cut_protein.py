from pymol import cmd, stored
import random
import warnings
import json
from glob import glob
from importlib import resources

def sort_residues(unsorted):
# sort residues and allow for non-integer resis
    intidx = []
    stridx = []
    for resi in unsorted:
        try:
            intidx.append(int(resi))
        except ValueError:
            stridx.append(resi)
    intidx.sort()

    insertmap = dict()
    for resi in stridx:
        resv = int(resi[:-1])
        letter = resi[-1]
        if resv not in insertmap.keys():
            insertmap[resv] = [letter]
        else:
            insertmap[resv].append(letter)


    for resv, letters in insertmap.items():
        letters.sort()
        idx = intidx.index(resv)
        insertion = ["%s%s" % (resv, letter) for letter in letters]
        intidx = intidx[:(idx + 1)] + insertion + intidx[(idx+1):]

    residues = [str(resi) for resi in intidx]
    return residues

# make dict of resis of each protein chain
def resisdict():
    resis = {}
    cmd.select("sidechains", "sidechain")
    stored.chains = cmd.get_chains("sidechains")
    for chain in stored.chains:
        stored.residues = []
        cmd.select("subsidechains", "sidechain and chain %s" % chain)
        cmd.iterate("subsidechains", "stored.residues.append(resi)")
        unsorted_residues = list(set(stored.residues))
        resis[f'{chain}'] = sort_residues(unsorted_residues)
    print('resis:',resis)
    return stored.chains, resis

def fragmentprotein(sub:str, monoC:str = None):
    """
    Fragments a protein between QM and MM regions

    Parameters
    ----------
    sub : str
        seed for growing QM region
        'ligand' or atom id number
        
    monoC : selection for the external charges
        Default: None

    NOTES

        - The ligand is monomer A, the QM region of the protein is
        monomer B, and the MM region of the protein is monomer C.

        - The QM region is initially cut such that peptide bonds are
        broken. The QM region then shrinks or expands such that only
        C_a - C_o bonds are broken.

        - Once a proper QM region is selected, we check for MM residues
        sandwiched between two QM regions. We move that residue from
        the MM region to the QM region.
    
    Returns
    -------
    None
    """
    name="fragmentprotein"

    out = open(glob('*.out')[0], 'a')

    # Get ligand info, add to system
    cmd.select("lig_A", "bm. hetatm and not sol. and not metals")
    
    # Carve away monomer C
    if monoC is not None:
        # Selecting by distance to ligand?
        if "be." in monoC:
            # Digest monoC declaration
            monoC = monoC.strip('"')
            cutoff = monoC.split()[-1]
            # Make pocket selection
            ## Select full residues or metal cofactors in the pocket
            # Cut between C (of carbonyl group) and N. 
            if sub == 'ligand':
                cmd.select('sys%s_B' % cutoff, "(byres not lig_A w. %s of lig_A) or (metals w. %s of lig_A) or (solvent w. %s of lig_A) " % (cutoff, cutoff, cutoff))
            else:
                sub = int(sub)
                cmd.select('sys%s_B' % cutoff, "(byres not lig_A w. %s of id %s) or (metals w. %s of id %s) or (solvent w. %s of id %s) " % (cutoff, sub, cutoff, sub, cutoff, sub))
            # Expand B such that only alpha carbon -- carbon bonds are broken. This only handles when N is on the QM side.
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem N) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem C) xt. 1 and elem O)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem Se) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem Se) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + (sys%s_B (xt. 1 and elem H))" % (cutoff, cutoff)) 
            ## Select everything else, that's monomer C
            cmd.select('mono_C', "not sys%s_B and not lig_A" % cutoff)
            # Expand C such that only alpha carbon -- carbon bonds are broken. This handles when N is on the MM side.
            cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 1 and elem S)")
            cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 2)")
            cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem H)")
            cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem N)")
            cmd.select("mono_C", "mono_C + ((mono_C and elem N) xt. 1)") 
            cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 1 and elem O)") 
            # Expand C so that endcaps aren't on fronteir regions
            cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn ACE)") 
            cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn NME)") 
            cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn NMA)") 
            # Now re-specify system B so that there are no overlapping atoms (in both system B and monoC)
            cmd.select("sys%s_B" % cutoff, "not mono_C and not lig_A")
            out.write('----------------------------------------------------------------------------------------------------\n')
            out.write('cut_protein'.center(100)+'\n')
            out.write('----------------------------------------------------------------------------------------------------\n')
            out.write('Making an initial cut:\n')
            out.write('Number of atoms in QM protein: ')
            out.write(f"{cmd.count_atoms('sys%s_B' % cutoff)}\n")
            # Identify any MM residues between two QM residues
            stored.Cresis = []
            # Get atoms in mono_C that are directly bound to system B
            cmd.select("boundary_cs", "mono_C and bound_to sys%s_B" % cutoff)
            # Add the residue numbers of these atoms to a list
            cmd.iterate("boundary_cs", "stored.Cresis.append(resi)")
            # For each residue bound directly to system B
            for r in stored.Cresis:
                resis = []
                up = str(int(r)+1)
                down = str(int(r)-1)
                # get atoms that are in the neighboring residue and part of system B and not C or O
                cmd.select("nextto", "sys%s_B and resi %s and not name C and not name O" % (cutoff, up))
                cmd.select("nextto", "nextto + (sys%s_B and resi %s and not name C and not name O)" % (cutoff, down))
                # find where there is more than one neighboring QM residue
                stored.neighbor_resis = []
                cmd.iterate("nextto", "stored.neighbor_resis.append(resi)")
                neighbor_resis = set(stored.neighbor_resis)
                if len(neighbor_resis) == 2:
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + resi %s" % (cutoff, r))
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem N) xt. 1)" % (cutoff, cutoff, r)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem C) xt. 1 and elem O)" % (cutoff, cutoff, r)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem S) xt. 1)" % (cutoff, cutoff, r)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + (br. (sys%s_B and elem S))" % (cutoff, cutoff)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem Se) xt. 1)" % (cutoff, cutoff, r)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and resi %s and elem Se) xt. 1)" % (cutoff, cutoff, r)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + (sys%s_B (xt. 1 and elem H))" % (cutoff, cutoff)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem N) xt. 1)" % (cutoff, cutoff)) 
                    cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem C) xt. 1 and elem O)" % (cutoff, cutoff)) 
                    out.write('Extending the QM region to avoid having an MM residue surrounded by QM residues on each side:\n')
                    out.write(f"Number of atoms in new QM region: {cmd.count_atoms('sys%s_B' % cutoff)}\n")
                    ## Select everything else, that's monomer C
                    cmd.select('mono_C', "not sys%s_B and not lig_A" % cutoff)
                    ## Exclude peptide bonds from border residues, i.e., B--C border is cutting across CA--(sidechain) bond
                    # Expand C such that only alpha carbon -- carbon bonds are broken. This handles when N is on the MM side.
                    cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 1 and elem S)")
                    cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 2)")
                    cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem H)")
                    cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem N)")
                    cmd.select("mono_C", "mono_C + ((mono_C and elem N) xt. 1)") 
                    cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 1 and elem O)") 
                    # Expand C so that endcaps aren't on fronteir regions
                    cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn ACE)") 
                    cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn NMA)") 
                    # Now re-specify system B so that there are no overlapping atoms (in both system B and monoC)
                    cmd.select("sys%s_B" % cutoff, "not mono_C and not lig_A")
                    cmd.hide("sticks", "mono_C")
                    cmd.show("sticks", "mono_C and sol.")
                    cmd.show("lines", "mono_C and not sol.")


    else:
        B_cutoff = ''
        # Make pocket selection
        cmd.select('sys%s_B' % B_cutoff, "not lig_A")
    
    cmd.save('QM.pdb', 'sys%s_B + lig_A' % cutoff)
    cmd.save('external.pdb', 'mono_C')
    cmd.save('ligand.pdb', 'lig_A')
    out.close()
    return

# make fragments of every 3 resn or every 1 resn
def makeresifragments(num_resis):
    # first, get ligand info, add to system
    cmd.select("lig_A", "bm. hetatm and not sol. and not metals")
    #charge_dict = make_charge_dict(f'{args.PDB_file}')
    stored.chains, resis = resisdict()
    #frag_charges = {}
    print('stored.chains', stored.chains)
    for chain in stored.chains:
        print(chain)
        resi_list = resis[chain]
        print('resi_list', resi_list)
        resi_len = len(resi_list)
        print('resi_len', resi_len)
        for i in range(resi_len):
            i = int(i)
            print('i',i)
            if num_resis == '3':
                fragment = resi_list[i:i+3]
                if len(fragment) == 3: #make sure every section has 3 residues (handles end cases)
                    cmd.select("frag", "chain %s and (resi %s or resi %s or resi %s)" % (chain, fragment[0], fragment[1], fragment[2]))
              #      charge = retrieve_charges(chain, fragment, charge_dict)
                    filename = f'{chain}_'+'-'.join(fragment)
              #      frag_charges[filename] = charge
                    cmd.save(f'{filename}.pdb',"frag",-1, 'pdb')
                    #frag_xyz = sele2xyz('frag', f'{chain}_'+'-'.join(fragment))
            elif num_resis == '1':
                fragment = [resi_list[int(i)]]
                if len(fragment) == 1: #make sure every section has 1 residue (handles end cases)
                    out = open(glob('*.out')[0], 'a')
                    cmd.select("sys%s_B" % fragment[0], "chain %s and (resi %s)" % (chain, fragment[0]))
                    # Expand B such that only alpha carbon -- carbon bonds are broken. This only handles when N is on the QM side.
                    cmd.select('sys%s_B' % fragment[0], "sys%s_B + ((sys%s_B and elem N) xt. 1)" % (fragment[0], fragment[0])) 
                    cmd.select('sys%s_B' % fragment[0], "sys%s_B + ((sys%s_B and elem C) xt. 1 and elem O)" % (fragment[0], fragment[0])) 
                    cmd.select('sys%s_B' % fragment[0], "sys%s_B + ((sys%s_B and elem Se) xt. 1)" % (fragment[0], fragment[0])) 
                    cmd.select('sys%s_B' % fragment[0], "sys%s_B + ((sys%s_B and elem Se) xt. 1)" % (fragment[0], fragment[0])) 
                    cmd.select('sys%s_B' % fragment[0], "sys%s_B + (sys%s_B (xt. 1 and elem H))" % (fragment[0], fragment[0])) 
                    # Expand B to include endcaps if next to QM region
                    cmd.select('sys%s_B' % fragment[0], "sys%s_B + ((sys%s_B and elem C) xt. 3 and resn ACE)" % (fragment[0], fragment[0]))
                    cmd.select('sys%s_B' % fragment[0], "sys%s_B + ((sys%s_B and elem C) xt. 3 and resn NMA)" % (fragment[0], fragment[0]))
                    cmd.select('sys%s_B' % fragment[0], "sys%s_B + ((sys%s_B and elem C) xt. 3 and resn NME)" % (fragment[0], fragment[0]))
                    #cmd.select('sys%s_B' % fragment[0], "sys%s_B + byres (sys%s_B and resn ACE)" % (fragment[0], fragment[0])) 
                    #cmd.select('sys%s_B' % fragment[0], "sys%s_B + byres (sys%s_B and resn NME)" % (fragment[0], fragment[0])) 
                    ## Select everything else, that's monomer C
                    cmd.select('mono_C', "not sys%s_B and not lig_A" % fragment[0])
                    # Expand C such that only alpha carbon -- carbon bonds are broken. This handles when N is on the MM side.
                    cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 1 and elem S)")
                    cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 2)")
                    cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem H)")
                    cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem N)")
                    cmd.select("mono_C", "mono_C + ((mono_C and elem N) xt. 1)") 
                    cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 1 and elem O)") 
                    # Expand C so that endcaps aren't on fronteir regions
                    #cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn ACE)") 
                    #cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn NMA)") 
                    # Now re-specify system B so that there are no overlapping atoms (in both system B and monoC)
                    cmd.select("sys%s_B" % fragment[0], "not mono_C and not lig_A")
                    out.write('----------------------------------------------------------------------------------------------------\n')
                    out.write('cut_protein'.center(100)+'\n')
                    out.write('----------------------------------------------------------------------------------------------------\n')
                    out.write('Making an initial cut:\n')
                    out.write('Number of atoms in QM protein: ')
                    out.write(f"{cmd.count_atoms('sys%s_B' % fragment[0])}\n")
                    filename = f'{chain}_{fragment[0]}'
                    cmd.save(f'{filename}.pdb','sys%s_B' % fragment[0],-1, 'pdb')
                    print('saved filename:', filename)
                    cmd.save(f'external_{filename}.pdb', 'mono_C')
                    cmd.save('ligand.pdb', 'lig_A')
                    out.close()
                    print('calling makepredict')
                    #cmd.do(f'makepredictionary {fragment[0]} {fragment[0]}')
                    print('entering makepredict explicit')
                    M1_atms = cmd.identify("(neighbor sys%s_B)" % fragment[0])
                    bond_dict = {}
                    bond_dict['QM'] = cmd.identify('sys%s_B' % fragment[0])
                    bond_dict['MM'] = cmd.identify('mono_C')
                    # making M1_x, Q1_x, M2_x, M3_x lists for each bond broken (x). 
                    for n,m in enumerate(M1_atms):
                        bond_dict[f'M1_{n+1}'] = [m]
                        bond_dict[f'Q1_{n+1}'] = cmd.identify(f'sys%s_B and neighbor id {m}' % fragment[0])
                        bond_dict[f'M2_{n+1}'] = cmd.identify(f'mono_C and neighbor id {m}')
                        bond_dict[f'M3_{n+1}'] = []
                        for x in bond_dict[f'M2_{n+1}']:
                            for a in cmd.identify(f'(mono_C and neighbor id {x}) and not id {m}'):
                                bond_dict[f'M3_{n+1}'].append(a)
                    # remove M1_x, M2_x, M3_x, and Q1_x atoms from MM and QM lists
                    M123_atms = []
                    for n, m in enumerate(M1_atms):
                        if f'M1_{n+1}' in bond_dict.keys():
                            for x in bond_dict[f'M1_{n+1}']:
                                M123_atms.append(x)
                        if f'M2_{n+1}' in bond_dict.keys():
                            for x in bond_dict[f'M2_{n+1}']:
                                M123_atms.append(x)
                        if f'M3_{n+1}' in bond_dict.keys():
                            for x in bond_dict[f'M3_{n+1}']:
                                M123_atms.append(x)
                        for x in bond_dict[f'Q1_{n+1}']:
                            if x in bond_dict['QM']:
                                bond_dict['QM'].remove(x)
                    bond_dict['MM'] = list(set(bond_dict['MM']).difference(M123_atms))
                    with open(f'pre-dictionary_{filename}.dat', 'w+') as dictfile:
                        dictfile.write(json.dumps(bond_dict))
        return

def makepredictionary(cutoff:str, filename:str='') -> None: 
    """ 
    Creates an initial version of the dictionary that assigns each atom
    to its region. Atoms in the fronteir region are named according to
    the bond cut as Q1_{bond}, M1_{bond}, M2_{bond}, M3_{bond}. All
    other QM atoms are added to QM, and all other MM atoms are added
    to MM, QM, MM, Q1_{bond}, M1_{bond}, etc.

    Parameters
    ----------
    cutoff: str
       radius in angstroms for making cut 

    Returns
    -------
    None
    """
    print('entering makepredict fxn')
    print('cutoff', cutoff)
    M1_atms = cmd.identify("(neighbor sys%s_B)" % cutoff)
    print('found M1_atms')
    bond_dict = {}
    bond_dict['QM'] = cmd.identify('sys%s_B' % cutoff)
    print('1:',bond_dict)
    bond_dict['MM'] = cmd.identify('mono_C')
    print('2:',bond_dict)
    # making M1_x, Q1_x, M2_x, M3_x lists for each bond broken (x). 
    for n,m in enumerate(M1_atms):
        bond_dict[f'M1_{n+1}'] = [m]
        print('3:',bond_dict)
        bond_dict[f'Q1_{n+1}'] = cmd.identify(f'sys%s_B and neighbor id {m}' % cutoff)
        print('4:',bond_dict)
        bond_dict[f'M2_{n+1}'] = cmd.identify(f'mono_C and neighbor id {m}')
        print('5:',bond_dict)
        bond_dict[f'M3_{n+1}'] = []
        for x in bond_dict[f'M2_{n+1}']:
            for a in cmd.identify(f'(mono_C and neighbor id {x}) and not id {m}'):
                bond_dict[f'M3_{n+1}'].append(a)
    # remove M1_x, M2_x, M3_x, and Q1_x atoms from MM and QM lists
    M123_atms = []
    for n, m in enumerate(M1_atms):
        if f'M1_{n+1}' in bond_dict.keys():
            for x in bond_dict[f'M1_{n+1}']:
                M123_atms.append(x)
        if f'M2_{n+1}' in bond_dict.keys():
            for x in bond_dict[f'M2_{n+1}']:
                M123_atms.append(x)
        if f'M3_{n+1}' in bond_dict.keys():
            for x in bond_dict[f'M3_{n+1}']:
                M123_atms.append(x)
        for x in bond_dict[f'Q1_{n+1}']:
            if x in bond_dict['QM']:
                bond_dict['QM'].remove(x)
    bond_dict['MM'] = list(set(bond_dict['MM']).difference(M123_atms))
    with open(f'pre-dictionary{filename}.dat', 'w+') as dictfile:
        dictfile.write(json.dumps(bond_dict))
    return

cmd.extend("fragmentprotein", fragmentprotein)
cmd.extend("resisdict", resisdict)
cmd.extend("makeresifragments", makeresifragments)
cmd.extend("makepredictionary", makepredictionary)

def run_cut_protein(pdb_file:str, sub:str, cutoff:str, num_resis:str = 0) -> None:
    """
   Calls the other necessary functions to cut the system in pdb_file
   by including everything that is {cutoff} angstroms away from {sub}

    Parameters
    ----------
    pdb_file: str
       path to pdb
    sub: str
       atom or group to grow QM region around
    cutoff: str
       radius in angstroms for making cut 

    Returns
    -------
    None
    """
    cmd.feedback("disable", "all", "everything")
    cmd.extend("makepredictionary", makepredictionary)
    cmd.reinitialize()
    cmd.load(pdb_file)
    with resources.path('sparcle_qc.data', 'cut_protein.py') as file_path:
        cut_path = str(file_path)
    cmd.do(f'run {cut_path}')
    if num_resis != '0':
        cmd.do(f'makeresifragments {num_resis}')
    else:
        cmd.do(f'fragmentprotein {sub}, monoC="be. {cutoff}"')
        cmd.do(f'makepredictionary {cutoff}')
