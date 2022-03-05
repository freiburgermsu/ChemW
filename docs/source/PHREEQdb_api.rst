PHREEQdb
++++++++++++++++++

The data environment, in a `Python IDE <https://www.simplilearn.com/tutorials/python-tutorial/python-ide>`_, is defined: 

.. code-block:: python

 import chemw
 phreeq_db = chemw.PHREEQdb(output_path = None, verbose = False, printing = False)

- *output_path* ``str``: optionally specifies an path to where the processed PHREEQ database file will be exported, where `None` selects the current working directory.
- *verbose* & *printing* ``bool``: optionally specifies whether progress or results of the calculations, respectively, are printed. The former is valuable for troubleshooting while the latter is beneficial for reviewing a readout summary of the calculations.

++++++++++
process()
++++++++++

A PHREEQ database file is processed into a JSON file of the elements and minerals, with their respective formula and MW: 

.. code-block:: python

 phreeq_db.process(db_path)

- *db_path* ``str``: The path to where the ``.dat`` PHREEQ database file that will be processed.


++++++++++++++++++++++++++
Accessible content
++++++++++++++++++++++++++
The ``PHREEQdb`` object retains numerous components that are accessible to the user: 

- *db_name* ``str``: The name of the database that is parsed in the ``process()`` function.
- *db*, *minerals*, & *elements* ``Pandas.DataFrame``: The entire PHREEQ database and the minerals and elements of the PHREEQ database, respectively, expressed in a Pandas Database object, and organized with labeled columns of the content. 
- *chem_mw* ``ChemMW``: An instance of the ``ChemMW`` object is loaded, which allows the user to access the ``ChemMW`` module through the ``PHREEQdb`` module.