from pymol import cmd
import random
import warnings
import json
from glob import glob

def pocketfrag(sub, split = "sidechains", monoC = None, prep_fisapt=False, cutcap=False):
    """
    DESCRIPTION

        Fragments a ligand@protein system and prepares F-/ISAPT input & accompanying fragment
        files.

    USAGE

        pocketfrag [, [split, [monoC, [prep_fisapt]]]]
    
    ARGUMENTS

        split : how to split the protein subsystem.
            Options: "residues", "sidechains" or "ss"
            Default: "sidechains"

        monoC : selection expression for defining ISAPT "monomer C", if any
            Default: None

        prep_fisapt : Boolean for whether or not to prepare F-/ISAPT input & fragment files
            Default: False

        cutcap : bool **EXPERIMENTAL** 
            Carve out the A:B:C subset from the entire protein by cutting bonds to anything
            else, and capping with the naive PyMOL cmd.h_add()
            Default: False

    NOTES

        * Specifics of the fragmentation scheme implemented herein are consistent with the
        original application of F-SAPT to protein-ligand binding by R. M. Parrish
        (Chem. Eur. J. 2017, 23, 7887–7890; doi: 10.1002/chem.201701031)
    
        * By default, the ligand is considered F-SAPT monomer A, the protein is considered
        F-SAPT monomer B, and any atoms satisfying the selection expression passed for the
        `monoC` argument are considered ISAPT monomer C.

        * A disulfide bridge between two cysteine residues is split into its own sidechain
        by cutting across the CA-//-SG bond, --S...S--, labeled according to the smaller
        residue ID (resi) of the two residues being bridged

        * If splitting on residues (`split="residues"`):
            * __**DANGER**__: Peptide bonds are cut between C and N, i.e., ...CA-C-//-N-CA...,
            which have partial double bond character!! This should only be done by advanced
            users who know _exactly_ what they're getting into!!

        * If splitting on sidechains:
            * atoms CA are included in the sidechain, not the backbone
            * peptide bonds are included as their own functional groups, labeled according
            to the smaller residue ID (resi) of the two residues being connected

        * If declaring monomer C using a "beyond `cutoff`" or "be. `cutoff`"  selection
        expression:
            * the full sidechain is included in the protein "pocket" for any residue with
            any atom within the cutoff distance of the ligand
            * peptide bonds connecting included sidechains and monomer C are _not_ included
            in the pocket, i.e., the single bond connecting the sidechain (monomer B) and
            backbone CA (monomer C) is cut.
    
    EXAMPLES

        # Fragment entire protein with respect to sidechains
        PyMOL> pocketfrag split="sidechain" 

        # Fragment protein wrt sidechains with an arbitrary selection expression defining
        # ISAPT monomer C
        PyMOL> pocketfrag split="ss", monoC="chain B and not sol."

        # Fragment protein wrt sidechains with cutoff distance from ligand defining ISAPT
        # monomer C
        PyMOL> pocketfrag split="ss", monoC="be. 10"

    """
    name="pocketfrag"

    # Warn about splitting over residues
    if split=='residues':
        warnings.warn("""DANGER!! Splitting over residues will cause peptide bonds with partial
double-bond character to be cut! This is _never_ advisable with F-/ISAPT
and may lead to garbage results in F-SAPT fragmentation!!""")

    out = open(glob('*.out')[0], 'a')


    # Get ligand info, add to system
    cmd.select("lig_A", "bm. hetatm and not sol. and not metals")
    
    # Now that ligand has been digested, color it
    # Carve away monomer C
    if monoC is not None:
        # Selecting by distance to ligand?
        if "be." in monoC and "w." not in monoC:
            # Digest monoC declaration
            monoC = monoC.strip('"')
            cutoff = monoC.split()[-1]
            # Make pocket selection
            ## Select full residues or metal cofactors in the pocket
            # CTS: Running command below cuts C (of carbonyl group) and N. Following lines adjust the cut
            if sub == 'ligand':
                cmd.select('sys%s_B' % cutoff, "(byres not lig_A w. %s of lig_A) or (metals w. %s of lig_A) or (solvent w. %s of lig_A) " % (cutoff, cutoff, cutoff))
            else:
                sub = int(sub)
                cmd.select('sys%s_B' % cutoff, "(byres not lig_A w. %s of id %s) or (metals w. %s of id %s) or (solvent w. %s of id %s) " % (cutoff, sub, cutoff, sub, cutoff, sub))
            # CTS expand B such that only alpha carbon -- carbon bonds are broken. This only handles when N is on the QM side.
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem N) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem C) xt. 1 and elem O)" % (cutoff, cutoff)) 
            #cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem S) xt. 1)" % (cutoff, cutoff)) 
            #cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem S) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem Se) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + ((sys%s_B and elem Se) xt. 1)" % (cutoff, cutoff)) 
            cmd.select('sys%s_B' % cutoff, "sys%s_B + (sys%s_B (xt. 1 and elem H))" % (cutoff, cutoff)) 
            ## Select everything else, that's monomer C
            cmd.select('mono_C', "not sys%s_B and not lig_A" % cutoff)
            ## Exclude peptide bonds from border residues, i.e., B--C border is cutting across CA--(sidechain) bond
            # CTS expand C such that only alpha carbon -- carbon bonds are broken. This handles when N is on the MM side.
            cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 1 and elem S)")
            cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 2)")
            cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem H)")
            cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem N)")
            cmd.select("mono_C", "mono_C + ((mono_C and elem N) xt. 1)") 
            cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 1 and elem O)") 
            # CTS: expand C so that endcaps aren't on fronteir regions
            cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn ACE)") 
            cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn NMA)") 
            # CTS now re-specify system B so that there are no overlapping atoms (in both system B and monoC)
            cmd.select("sys%s_B" % cutoff, "not mono_C and not lig_A")
            out.write('----------------------------------------------------------------------------------------------------\n')
            out.write('cut_protein'.center(100)+'\n')
            out.write('----------------------------------------------------------------------------------------------------\n')
            out.write('Making an initial cut:\n')
            out.write('Number of atoms in QM protein: ')
            out.write(f"{cmd.count_atoms('sys%s_B' % cutoff)}\n")
            cmd.hide("sticks", "mono_C")
            cmd.show("sticks", "mono_C and sol.")
            cmd.show("lines", "mono_C and not sol.")
            # CTS identify any MM residues between two QM residues
            stored.Cresis = []
            # Get atoms in mono_C that are directly bound to system B
            cmd.select("boundary_cs", "mono_C and bound_to sys%s_B" % cutoff)
            # Add the residue numbers of these atoms to a list
            cmd.iterate("boundary_cs", "stored.Cresis.append(resi)")
            #print(stored.Cresis)
            # For each residue bound directly to system B
            for r in stored.Cresis:
                resis = []
                up = str(int(r)+1)
                down = str(int(r)-1)
                # get atoms that are in the neighboring residue and part of system B and not C or O
                cmd.select("nextto", "sys%s_B and resi %s and not name C and not name O" % (cutoff, up))
                cmd.select("nextto", "nextto + (sys%s_B and resi %s and not name C and not name O)" % (cutoff, down))
                #cmd.save('%s.pdb' % r, "nextto")
                # find where there is more than one neighboring QM residue
                stored.neighbor_resis = []
                cmd.iterate("nextto", "stored.neighbor_resis.append(resi)")
                neighbor_resis = set(stored.neighbor_resis)
                if len(neighbor_resis) == 2:
                    #print(neighbor_resis)
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
                    # CTS expand C such that only alpha carbon -- carbon bonds are broken. This handles when N is on the MM side.
                    cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 1 and elem S)")
                    cmd.select("mono_C", "mono_C + ((mono_C and elem S) xt. 2)")
                    cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem H)")
                    cmd.select("mono_C", "mono_C + ((mono_C and name CA) xt. 1 and elem N)")
                    cmd.select("mono_C", "mono_C + ((mono_C and elem N) xt. 1)") 
                    cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 1 and elem O)") 
                    # CTS: expand C so that endcaps aren't on fronteir regions
                    cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn ACE)") 
                    cmd.select("mono_C", "mono_C + ((mono_C and elem C) xt. 3 and resn NMA)") 
                    # CTS now re-specify system B so that there are no overlapping atoms (in both system B and monoC)
                    cmd.select("sys%s_B" % cutoff, "not mono_C and not lig_A")
                    cmd.hide("sticks", "mono_C")
                    cmd.show("sticks", "mono_C and sol.")
                    cmd.show("lines", "mono_C and not sol.")
            # loop through each residue at boundary (in monoC)
            # expand; if expansion returns two atoms in QM
        elif "be." in monoC and "w." in monoC:
            monoC = monoC.strip('"').split()
            B_cutoff = monoC[monoC.index('be.') + 1]
            C_cutoff = monoC[monoC.index('w.') + 1]
            multilevel(B_cutoff, C_cutoff, cutcap=bool(cutcap))
        else:
            cmd.select('mono_C', monoC.strip('"'))


    else:
        B_cutoff = ''
        # Make pocket selection
        cmd.select('sys%s_B' % B_cutoff, "not lig_A")
    
    cmd.save('QM.pdb', 'sys%s_B + lig_A' % cutoff)
    cmd.save('external.pdb', 'mono_C')
    cmd.save('ligand.pdb', 'lig_A')
    out.close()
    return

def make_dictionary(cutoff): #for capping with Caroline's code
    M1_atms = cmd.identify("(neighbor sys%s_B)" % cutoff)
    bond_dict = {}
    bond_dict['QM'] = cmd.identify('sys%s_B' % cutoff)
    bond_dict['MM'] = cmd.identify('mono_C')
    # making M1_x, Q1_x, M2_x, M3_x lists for each bond broken (x). 
    for n,m in enumerate(M1_atms):
        bond_dict[f'M1_{n+1}'] = [m]
        bond_dict[f'Q1_{n+1}'] = cmd.identify(f'sys%s_B and neighbor id {m}' % cutoff)
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
    with open('pre-dictionary.dat', 'w+') as dictfile:
        dictfile.write(json.dumps(bond_dict))
    return

cmd.extend("pocketfrag", pocketfrag)
cmd.extend("make_dictionary", make_dictionary)

def run_cut_protein(pdb_file, sub, cutoff):
    cmd.feedback("disable", "all", "everything")
    cmd.reinitialize()
    cmd.load(pdb_file)
    cmd.do('run ../cut_protein.py')
    cmd.do(f'pocketfrag {sub}, split="ss", monoC="be. {cutoff}"')
    cmd.do(f'make_dictionary {cutoff}')
