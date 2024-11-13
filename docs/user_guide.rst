User Guide
===============

General Information
*********************

.. dropdown:: Cutting and Capping Bonds
    
    Alpha carbon - carbonyl carbon bonds are cut to separate the QM region from the MM region. The valencies on the QM atoms of the cut bonds are satisfied by placing a hydrogen along the cut bond.

.. dropdown:: QM/MM Boundary Charges

    Point charges too close to the capped QM region may cause overpolarization. Users have the option of choosing one of nine charge schemes to alter charges at this boundary. These schemes are outlined below:

.. dropdown:: Resulting Input File for QM Calculation

    The result of running Sparcle_QC is an input file for either Psi4, Q-Chem, or NWChem. In general, the input file will look similar to the files below.

    .. tab-set::
    
        .. tab-item:: Psi4
    
            .. code-block:: 
    
                int main(const int argc, const char **argv) {
                    return 0;
                }
    
        .. tab-item:: Q-Chem
    
            .. code-block:: 
    
                def main():
                    return
    
        .. tab-item:: NWChem
    
            .. code-block:: 
    
                class Main {
                    public static void main(String[] args) {
                    }
                }

Options
********

.. dropdown:: Required
    
    .. csv-table:: 
        :file: options/required.csv

.. dropdown:: Force Fields

    .. csv-table:: 
        :file: options/forcefields.csv

.. dropdown:: Psi4

    .. csv-table:: 
        :file: options/psi4.csv

.. dropdown:: Q-Chem

    .. csv-table:: 
        :file: options/qchem.csv

.. dropdown:: NWChem

    .. csv-table:: 
        :file: options/nwchem.csv

.. dropdown:: Other

    .. csv-table:: 
        :file: options/other.csv


Example Inputs
***************

.. dropdown:: F-SAPT with Psi4 and Amber

    The following input file will create an F-SAPT file for the protein:ligand complex to be run with Psi4. It will also create the functional group partitions needed for post-processing, fA.dat and fB.dat. 
    
    .. code-block::
 
        pdb_file: 2cji.pdb
        pre-capped: true
        cutoff: 5
        seed: ligand
        charge_scheme: BRC
        ligand_charge: 0
        method: fisapt0
        fisapt_partition: true
        basis_set: aug-cc-pv(D+d)z
        amber_ff: ff19SB
        env_path: /user/miniconda3/envs/sparcle_qc/
        water_model: opc
        o_charge: 0
        h_charge: 0.6791
        ep_charge: -1.3582
        software: psi4
        mem: 60 GB
        nthreads: 10


.. dropdown:: B3LYP with Q-Chem and CHARMM

    The following input file will create prepare 3 Q-Chem files: one with the ligand  (fully QM), one with the protein (QM/MM), and one with the complex (QM/MM). These could be used to calculate a supermolecular interaction energy. We will turn on counterpoise correction, which will include ghost atoms for the QM dimer in all three files.
    
    .. code-block::
 
        pdb_file: 3qxp.pdb
        cutoff: 5
        seed: ligand
        charge_scheme: DZ3
        ligand_charge: 0
        method: b3lyp
        basis_set: 6-31G*
        charmm_rtf: top_all36_prot.rtf
        charmm_prm: par_all36m_prot.prm
        water_model: tip3p
        software: q-chem

.. dropdown:: HF with NWChem and Amber

    The following input file will create prepare 3 NWChem files: one with the ligand  (fully QM), one with the protein (QM/MM), and one with the complex (QM/MM). These could be used to calculate a supermolecular interaction energy. We will turn on counterpoise correction, which will include ghost atoms for the QM dimer in all three files. The QM region will grow starting from a single ligand atom.
    
    .. code-block::

        pdb_file: 2cji.pdb
        pre-capped: true
        cutoff: 8.5
        seed: 4247
        seed_file: 4yff.pdb
        charge_scheme: BRC
        ligand_charge: 0
        method: hf 
        basis_set: aug-cc-pv(D+d)z
        amber_ff: ff19SB
        env_path: /usr/miniconda3/envs/emb_sapt/
        water_model: opc
        o_charge: 0
        h_charge: 0.6791
        ep_charge: -1.3582
        software: nwchem
        nwchem_scratch: /scratch/user
        nwchem_perm: /scratch/user
        mem: 60 GB

.. dropdown:: Templating a QM Region for Congeneric Ligands

    Studies that compare a protein with two similar ligand structures may choose to equilibrate protein structures for each ligand. In this case, the two PDBs may be similar in structure, but not identical, and their coordinates likely will not match. Here, we show the steps of (1) creating a SAPT input file for one ligand (named methyl), then (2) using the QM region of 'methyl' as a template for cutting the QM region of the other ligand, named 'chlorine'.


    Step 1, the following is 'methyl.in':

    .. code-block::

        pdb_file: 2cji_methyl.pdb
        pre-capped: true
        cutoff: 5
        seed: ligand
        charge_scheme: BRC
        ligand_charge: 0
        method: fisapt0
        basis_set: aug-cc-pv(D+d)z
        amber_ff: ff19SB
        env_path: /user/miniconda3/envs/sparcle_qc/
        water_model: opc
        o_charge: 0
        h_charge: 0.6791
        ep_charge: -1.3582
        software: psi4
        mem: 60 GB
        nthreads: 10
 
       
    Step 2, the following is 'chlorine.in':

    .. code-block::

        pdb_file: 2cji_chlorine.pdb
        pre-capped: true
        template_path: methyl/cx_autocap_fixed.pdb
        charge_scheme: BRC
        ligand_charge: 0
        method: fisapt0
        basis_set: aug-cc-pv(D+d)z
        amber_ff: ff19SB
        env_path: /user/miniconda3/envs/sparcle_qc/
        water_model: opc
        o_charge: 0
        h_charge: 0.6791
        ep_charge: -1.3582
        software: psi4
        mem: 60 GB
        nthreads: 10
  
    With the two SAPT files, a relative interaction energy can be computed, giving insight into which ligand is more stable within the protein pocket. 
      
.. dropdown:: Convergence Study with Increasing QM Region Size via API

    A Python loop can be used to generate multiple input files with an increasing size of the QM region. We can increase the size of the QM region by incrementing the cutoff. An example Python script is below.    

    .. code-block::
 
        import sparcle_qc

	inputs = {
            'pdb_file': '2cji.pdb',
            'pre-capped': 'True',
            'seed': 'ligand',
            'charge_scheme': 'BRC',
            'ligand_charge': 0,
            'method': 'fisapt0',
            'basis_set': 'aug-cc-pv(D+d)z',
            'amber_ff': 'ff19SB',
            'env_path': '/usr/miniconda3/envs/emb_sapt/',
            'water_model': 'opc' ,
            'o_charge': 0,
            'h_charge': 0.6791,
            'ep_charge': -1.3582,
            'software': 'psi4',
            'mem': '60 GB',
            'nthreads': 10}

	cutoffs = [3, 4, 5]

	for c in cutoffs:
	    inputs['cutoff'] = f'{c}'
	    inputs['input_filename'] = f'cutoff_{c}.in'
	    print(inputs)
	    sparcle_qc.run_sparcle(user_options = inputs)

                                                     
