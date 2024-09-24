from pymol.cgo import *
from pymol import cmd, editor


def fix_numbers_amber(pdb_file: str) -> None:
    """
    When given a pdb, creates a new copy {pdb_file}_fixed.pdb that has the protein residues followed by waters and then the ligand
    Corrects for any mistakes in atom or residue numbering that may have been caused by manipulation of the system in pymol
    Ensures that the ligand atoms are labeled as HETATM

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
            out.write(f'ATOM  {atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}           {line[66:len(line)].strip():<3}\n')
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
            out.write(f'ATOM  {atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{line[66:len(line)].strip():>12}\n')
    for line in ligand_lines:
        if len(line)>70 and line[0:6].strip()=='ATOM' or line[0:6].strip()=='HETATM':
            atomnum +=1
            if line[22:26].strip()!=oldres:
                resnum+=1
                oldres = line[22:26].strip()
            out.write(f'HETATM{atomnum:>5}{line[11:16].strip():>5}{line[16:20].strip():>4}{line[20:22].strip():>2}{resnum:>4}{line[30:38].strip():>12}{line[38:46].strip():>8}{line[46:54].strip():>8}{line[54:60].strip():>6}{line[60:66].strip():>6}{line[66:len(line)].strip():>12}\n')
    if 'cx' in pdb_file:
        out.write('CONECT\n')
    out.write('END')
    out.close()
def autocap(pdb_file: str) -> None:
    """
    When given an amber complex pdb, caps the ends of each chain with a capping group (NME on C-terminal and ACE on N-terminal)
    and then renames the default pymol atom names to amber atom names for the new capping residues
    Saves the entire new, capped complex as cx_autocap.pdb, just the capped protein as prot_autocap.pdb, and the ligand as ligand.pdb

    Parameters
    ----------
    pdb_file: str
        path to uncapped complex pdb

    Returns
    -------
    None
    """
    cmd.reinitialize()
    
    # Load PDB
    cmd.load(pdb_file,"pdb")
    cmd.show("sticks", "all")
    cmd.label("all", "name")
    
    cmd.select("sidechains", "sidechain")
    chains = []
    chains = cmd.get_chains("sidechains")
    for chain in chains:
        cmd.select('capped_n', 'not name H1 and element H and bound_to (first name N and chain %s)' % chain)
        cmd.remove('capped_n')
        cmd.select("capped_c", "last name OXT and chain %s" % chain)
        cmd.remove("capped_c")
        cmd.alter('first name H1 and chain %s' % chain, "name='H'")
        #cmd.set('retain_order', 0)
        editor.attach_amino_acid("last name C and chain %s" % chain, 'nme' )
        editor.attach_amino_acid("first name N and chain %s" % chain, 'ace')
    
    
    
    cmd.alter('resname NME and name 1HH3', "name='H1'")
    cmd.alter('resname NME and name 2HH3', "name='H2'")
    cmd.alter('resname NME and name 3HH3', "name='H3'")
    cmd.alter('resname NME and name CH3', "name='C'")
    
    
    cmd.alter('resname ACE and name 1HH3', "name='H1'")
    cmd.alter('resname ACE and name 2HH3', "name='H2'")
    cmd.alter('resname ACE and name 3HH3', "name='H3'")
    
    cmd.save(f"cx_autocap.pdb", "pdb")
    
    cmd.select("ligand", "all and not polymer and not metals and not solvent and not resn nme and not resn ace")
    cmd.save("ligand.pdb", "ligand")
    cmd.remove("ligand")
    
    cmd.save(f"prot_autocap.pdb", "pdb")
	# creating a cpptraj script for the given pdb file
def skip_autocap(pdb_file: str) -> None:
    """
    When given an amber complex pdb that already has capping residues at each termini, this function is called instead of autocap(pdb_file) 
    Saves just the capped protein as prot_autocap.pdb and the ligand as ligand.pdb

    Parameters
    ----------
    pdb_file: str
        path to capped complex pdb

    Returns
    -------
    None
    """
    cmd.reinitialize()
    # Load PDB
    
    cmd.load(pdb_file,"pdb")
    cmd.show("sticks", "all")
    cmd.label("all", "name")
    
    cmd.select("ligand", "all and not polymer and not metals and not solvent and not resn nme and not resn ace")
    cmd.save("ligand.pdb", "ligand")
    cmd.remove("ligand")
    
    cmd.save(f"prot_autocap.pdb", "pdb")
def write_cpptraj(pdb_file: str)-> None:
    """
    Writes a cpptraj input file for the provided pdb

    Parameters
    ----------
    pdb_file: str
        path to complex pdb

    Returns
    -------
    None
    """
    with open('cpptraj.in', 'w') as f:
        f.write(f'parm {pdb_file}\n')
        f.write(f'loadcrd {pdb_file} name tmp1\n')
        f.write(f'prepareforleap crdset tmp1 name tmp2 pdbout uncapped.pdb nosugar\n')
def write_cpptraj_skip_autocap(pdb_file: str) -> None:
    """
    Writes a cpptraj input file for the provided pdb, but with a different name than write_cpptraj since the complex is already capped

    Parameters
    ----------
    pdb_file: str
        path to complex pdb

    Returns
    -------
    None
    """
    with open('cpptraj.in', 'w') as f:
        f.write(f'parm {pdb_file}\n')
        f.write(f'loadcrd {pdb_file} name tmp1\n')
        f.write(f'prepareforleap crdset tmp1 name tmp2 pdbout cx_autocap.pdb nosugar\n')
def write_tleap(forcefield: str, water_model: str) -> None:
    """
    Writes a tleap input file that loads the given forcefield and water model

    Parameters
    ----------
    forcefield: str
        name of amber forcefield (ex. ff19SB)
    water_model: str
        name of the desired water model (ex. OPC)

    Returns
    -------
    None
    """
    with open('tleap.in', 'w') as f:
        f.write(f'source leaprc.protein.{forcefield}\n')
        f.write(f'source leaprc.water.{water_model}\n')
        f.write(f'mol = loadPdb "prot_autocap_fixed.pdb"\n')
        f.write(f'savemol2 mol prot_autocap_fixed.mol2 1\n')
        f.write('quit')
