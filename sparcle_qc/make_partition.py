import MDAnalysis as mda

def partition(pdb_file: str) -> None:
    """
    Makes a partition file for each monomer for use in fisapt in psi4 (fA.dat and fB.dat)
    All ligand atoms are in a single functional group 
    Protein side chains, waters, and ions are in a functional group labeled by the resname+resnum (eg. ALA250)
    Peptide bond atoms from neighboring residues (C,O,N,H) are in a functional group labeled by PEP+lower_residue_number 
    (eg. the C and O from residue 250 and N and H from residue 251 would be labeled PEP250)

    Parameters
    ----------
    pdb_file: str
        path to complex pdb

    Returns
    -------
    None
    """
    #this function takes in a pdb file and creates the necessary partition files for FSAPT
    #fA.dat is the partition for the ligand
    #fB.dat is the partition of the protein and waters
    u = mda.Universe(pdb_file)
    fA = {}
    fB = {}
    # CHARMM terminal patch residues from https://charmm-gui.org/?doc=help&id=terminal_patch
    prot = u.select_atoms("((protein or resname ACE or resname NME or resname NTER or resname NNEU or resname GLYP or resname PROP or resname ACP or resname CTER or resname CNEU or resname CT1 or resname CT2 or resname CT3 or resname 5TER or resname 3TER) and not name C and not name O and not name N and not name H and not resname CYX and not resname CYS) or resname LIN or resname HOH or resname TIP or resname TIP3")
    last_prot_atom = int(prot[-1].id)
    for atom in prot:
        if f"{atom.resname}{atom.resnum}" not in fA.keys():
            fA[f"{atom.resname}{atom.resnum}"] = str(atom.id)
        else:
            fA[f"{atom.resname}{atom.resnum}"] = fA[f"{atom.resname}{atom.resnum}"] + " " +  str(atom.id)

    #Detecting which atoms are in disulfide bonds the two bridged SG atoms in their own functional group (DIS) 
    sg_atoms = u.select_atoms("protein and name SG")
    for atom in sg_atoms:
        sg_num = atom.id
        print('sg_num:', sg_num)
        other_sg = u.select_atoms(f'protein and name SG and around 3 id {sg_num}')
        close_h = u.select_atoms(f'protein and (name HG or name HG1) and around 1.5 id {sg_num}')
        print('other_sg:', other_sg)
        if len(other_sg) > 0 and len(close_h) == 0:
            cur_resnum = atom.resid
            print('cur_resnum:', cur_resnum)
            other_resnum = other_sg.resids[0]
            if f"DIS{other_resnum}" not in fA.keys():
                fA[f"DIS{cur_resnum}"] = str(atom.id)
                fA[f"CYX{cur_resnum}"] = ''
                for atom in u.select_atoms(f'resnum {cur_resnum}'):
                    if atom.name != 'SG':
                        fA[f"CYX{cur_resnum}"] += " " + str(atom.id)
            else:
                fA[f"DIS{other_resnum}"] = fA[f"DIS{other_resnum}"] + " " +  str(atom.id)
                fA[f"CYX{other_resnum}"] = ''
                for atom in u.select_atoms(f'resnum {other_resnum}'):
                    if atom.name != 'SG':
                        fA[f"CYX{other_resnum}"] += " " + str(atom.id)
        else: #then not a disulfide bridge
            cur_resnum = atom.resid
            print('else cur_resnum:', cur_resnum)
            fA[f"CYS{cur_resnum}"] = str(atom.id)
            for atom in u.select_atoms(f'resnum {cur_resnum}'):
                fA[f"CYS{cur_resnum}"] += " " + str(atom.id)
    #creating pep functional groups containing atoms of two neighboring residues 
    pep = u.select_atoms("protein or resname ACE or resname NME or resname NTER or resname NNEU or resname GLYP or resname P    ROP or resname ACP or resname CTER or resname CNEU or resname CT1 or resname CT2 or resname CT3 or resname 5TER or resname 3TER")
    pep = pep.select_atoms("name C or name O or name N or name H or name HN")
    for atom in pep:
        if atom.name =="N" or atom.name == "H" or atom.name == "HN":
            pepnum = int(atom.resnum)-1
            if f"PEP{pepnum}" not in fA.keys():
                fA[f"PEP{pepnum}"] = str(atom.id)
            else:
                fA[f"PEP{pepnum}"] = fA[f"PEP{pepnum}"] + " " +  str(atom.id)
        else:

            if f"PEP{atom.resnum}" not in fA.keys():
                fA[f"PEP{atom.resnum}"] = str(atom.id)
            else:
                fA[f"PEP{atom.resnum}"] = fA[f"PEP{atom.resnum}"] + " " +  str(atom.id)
    
    #selecting ligand and making dictionary
    ligand = u.select_atoms('not protein and not resname HOH and not resname TIP and not resname TIP3 and not resname LIN and not resname ACE and not resname NME and not resname NTER and not resname NNEU and not resname GLYP and not resname PROP and not resname ACP and not resname CTER and not resname CNEU and not resname CT1 and not resname CT2 and not resname CT3 and not resname 5TER and not resname 3TER')

    ligand_count = 0
    for atom in ligand:
        ligand_count +=1
        atom_num = ligand_count+last_prot_atom

        if f"LIGAND" not in fB.keys():
            fB[f"LIGAND"] = str(atom.id)
        else:
            fB[f"LIGAND"] = fB[f"LIGAND"] + " " +  str(atom.id)
    
    #writing results to partition files
    with open('fB.dat', "x") as f:
        for residue in fA:
            f.write(f"{residue} {fA[residue]}\n")
    with open("fA.dat","x") as f:
        for ligand in fB:
            f.write(f"{ligand} {fB[ligand]}\n")

    return


