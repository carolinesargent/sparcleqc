Getting Started
===============
Installing Sparcle_QC
---------------------
To get started, you will need to install sparcle_qc and its dependencies. We have provided a yaml file to create a conda environment containing all of the necessary packages.

First, clone the `Sparcle_QC repository <https://github.com/carolinesargent/sparcle_qc>`_

Change your working directory to sparcle_qc/:

.. code-block:: bash

    cd sparcle_qc

Create a conda environment from the yaml file provided in this repository:

.. code-block:: bash

    conda env create -f sparcle_qc.yaml

Command Line Usage 
------------------
Sparcle_QC is now installed in this enviornment and can be called on the command line in any directory using the following syntax:

.. code-block:: bash

    sparcle_qc input_file.in

Python API
----------    
Alternatively, for more advanced scripting utilities, Sparcle_QC can be imported as a python package:

.. code-block:: python
    
    import sparcle_qc

Sparcle_QC can be called by referencing an input file:

.. code-block:: python

    sparcle_qc.run_sparcle(input_file = 'input_file.in')

or by passing a dictionary of inputs: 

.. code-block:: python

    inputs = {
    'input_filename': 'test14.in',
    'pdb_file': '3QXP_templated_amber.pdb',
    'cutoff': 6,
    'seed': 'ligand',
    'charge_scheme': 'Z2',
    'ligand_charge': 0,
    'method': 'hf',
    'basis_set': 'aug-cc-pv(D+d)z',
    'amber_ff': 'ff19SB',
    'env_path': f'{os.getcwd()}/',
    'water_model': 'opc' ,
    'o_charge': 0,
    'h_charge': 0.6791,
    'ep_charge': -1.3582,
    'software': 'psi4',
    'mem': '100 GB',
    'nthreads': 16,
    }

    sparcle_qc.run_sparcle(user_options = inputs)


Examples 
--------

For more details on the inputs to Sparcle_QC and example inputs, check out the :doc:`user_guide`.

