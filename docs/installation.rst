Installing SparcleQC
=====================
Creating a SparcleQC Enviornment
---------------------------------
To get started, you will need to install sparcle_qc and its dependencies. We have provided a yaml file to create a conda environment containing all of the necessary packages.

First, clone the `SparcleQC repository <https://github.com/carolinesargent/sparcle_qc>`_:

.. code-block:: bash

   git clone https://github.com/carolinesargent/sparcle_qc.git

Change your working directory to sparcle_qc/:

.. code-block:: bash

    cd sparcle_qc

Create a conda environment from the yaml file provided in this repository:

.. code-block:: bash

    conda env create -f sparcle_qc.yaml

Command Line Usage 
------------------
SparcleQC is now installed in this enviornment and can be called on the command line in any directory using the following syntax:

.. code-block:: bash

    sparcle_qc input_file.in

Python API
----------    
Alternatively, for more advanced scripting utilities, SparcleQC can be imported as a python package:

.. code-block:: python
    
    import sparcle_qc

SparcleQC can be called by referencing an input file:

.. code-block:: python

    sparcle_qc.run_sparcle(input_file = 'input_file.in')

or by passing a dictionary of inputs: 

.. code-block:: python

    inputs = {
    'input_filename': 'example.in',
    'pdb_file': 'complex_pdb_file.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'Z2',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'software': 'psi4',
    }

    sparcle_qc.run_sparcle(user_options = inputs)


Running SparcleQC 
------------------

To learn more about SparcleQC, check out :doc:`getting_started`.

