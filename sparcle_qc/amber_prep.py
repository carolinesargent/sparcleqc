from pymol.cgo import *
from pymol import cmd, editor
import argparse

def autocap(pdb_file):
    
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
def write_cpptraj(pdb_file):
    with open('cpptraj.in', 'w') as f:
        f.write(f'parm {pdb_file}\n')
        f.write(f'loadcrd {pdb_file} name tmp1\n')
        f.write(f'prepareforleap crdset tmp1 name tmp2 pdbout uncapped.pdb nosugar\n')
def write_cpptraj_skip_autocap(pdb_file):
    with open('cpptraj.in', 'w') as f:
        f.write(f'parm {pdb_file}\n')
        f.write(f'loadcrd {pdb_file} name tmp1\n')
        f.write(f'prepareforleap crdset tmp1 name tmp2 pdbout cx_autocap.pdb nosugar\n')
def write_tleap(forcefield, water_model):
    with open('tleap.in', 'w') as f:
        f.write(f'source leaprc.protein.{forcefield}\n')
        f.write(f'source leaprc.water.{water_model}\n')
        f.write(f'mol = loadPdb "prot_autocap_fixed.pdb"\n')
        f.write(f'savemol2 mol prot_autocap_fixed.mol2 1\n')
        f.write('quit')
def skip_autocap(pdb_file):
    cmd.reinitialize()
    # Load PDB
    
    cmd.load(pdb_file,"pdb")
    cmd.show("sticks", "all")
    cmd.label("all", "name")
    
    cmd.select("ligand", "all and not polymer and not metals and not solvent and not resn nme and not resn ace")
    cmd.save("ligand.pdb", "ligand")
    cmd.remove("ligand")
    
    cmd.save(f"prot_autocap.pdb", "pdb")
