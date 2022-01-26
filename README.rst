Calculating the Molecular Weight of a Chemical
____________________________________________________________________

|PyPI version| |Actions Status| |Downloads| |License|

.. |PyPI version| image:: https://img.shields.io/pypi/v/chemw.svg?logo=PyPI&logoColor=brightgreen
   :target: https://pypi.org/project/chemw/
   :alt: PyPI version

.. |Actions Status| image:: https://github.com/freiburgermsu/chemw/workflows/Test%20ChemW/badge.svg
   :target: https://github.com/freiburgermsu/chemw/actions
   :alt: Actions Status

.. |License| image:: https://img.shields.io/badge/License-MIT-blue.svg
   :target: https://opensource.org/licenses/MIT
   :alt: License

.. |Downloads| image:: https://pepy.tech/badge/chemw
   :target: https://pepy.tech/project/chemw
   :alt: Downloads


The molecular weight (MW) can be algebraically calculated from any chemical formula that adheres to `chemical conventions <https://en.wikipedia.org/wiki/Chemical_formula>`_, which is exemplified in subsequent examples. 

The `ChemMW` object of ``ChemW`` parses a chemical formula string -- which consists of any combination of elements and decimal stoichiometry -- and precisely calculates the MW of the chemical formula, based upon the current physical accuracy that is embedded in the ``periodic_table`` dictionary of the ``chemicals`` module.

The `PHREEQdb` object of ``ChemW`` parses a `PHREEQ database <https://www.usgs.gov/software/phreeqc-version-3>`_ via the `ChemMW` object and exports a JSON of mineral masses for all of the described minerals in the database. This unique application of the `ChemMW` object has been applied as the pivotal means of predicting the mass of mineral scaling in the `ROSSpy module <https://pypi.org/project/ROSSpy/>`_ for reverse osmosis research.

The ``ChemW`` module is offered with the `MIT License <https://opensource.org/licenses/MIT>`_\. Examples of the module are available in the examples directory of the `ChemW GitHub repository <https://github.com/freiburgermsu/ChemW>`_. Please submit errors or inaccuracies as `GitHub issues <https://github.com/freiburgermsu/ChemW/issues>`_ so that they may be resolved.

+++++++++++++
installation
+++++++++++++

The following command installs ``ChemW`` in a command prompt/terminal environment::
 
 pip install chemw

_________________

ChemMW
++++++++++++++++++

+++++++++++
__init__
+++++++++++

The data environment, in a `Python IDE <https://www.simplilearn.com/tutorials/python-tutorial/python-ide>`_, is defined: 

.. code-block:: python

 import chemw
 chem_mw = chemw.ChemMW(verbose = False, printing = True)

- *verbose* & *printing* ``bool``: specifies whether troubleshooting information or MW results will be printed, respectively.

++++++++++++++++
mass()
++++++++++++++++

The parameterized data is fitted to the Hill equation, with the following arguments and their default values:

.. code-block:: python

 chem_mw.mass(formula)

- *formula* ``str``: parameterizes the chemical formula for which the MW is desired. The acceptable formats for the formula are quite broad, which are exemplified in the following formulae:

===================================================  ===================================================================================
 Example chemical                                      Format option
===================================================  ===================================================================================
 ``'C60_H120_O2'``                                       Underscores can arbitrarily separate content, since these are ignored by ``ChemMW``.
``'Na2.43_Cl_(OH)2_(OH)1.2_(OH)'``                      An arbitrary number of groups can be distinguished in the chemical formula, 
                                                            with ``()`` denoting the boundaries of the specified group.
  ``'Na2.43Cl(Ca(OH)2)1.2'``                             Chemical groups can be nested, with differing stoichiometric values.
 ``'Na2.43Cl:2H2O'``                                     Water molecules can be complexed, 
                                                               with a leading stoichiometric quantity of the complexation.
``'Na.96Al.96Si2.04O6:H2O'``                            Stoichiometry can be any decimal for any atom in a molecule, 
                                                                and even omit a leading zero.
``'Na2SO4:3K2SO4'``                                              Non-water entities can be complexed.
``'CaCl2:(MgCl2)2:12H2O'``                              Multiple complexations can be applied with repeated ``:`` separators. 
 ``'Ca1.019Na.136K.006Al2.18Si6.82O18:7.33H2O'``       The complexity, while remaining within the aforementioned format, is arbitrary.
===================================================  ===================================================================================


++++++++++++++++++++++++++
Accessible content
++++++++++++++++++++++++++
The ``ChemMW`` object retains numerous components that are accessible to the user: 

- *mw* ``float``: The MW of the parameterized chemical formula.
- *proportions* ``dict``: The ratio of elements in the chemical formula. This loses accuracy with the grouped elements, and is being improved.
- *formula* ``str``: The original chemical formula as a string.
- *groups* ``int``: A numerical counter for the quantity of chemical groups that are 
- *group_masses* ``dict``: A dictionary for the masses of each nesting level in a molecule.


_________________

PHREEQdb
++++++++++++++++++


++++++++++
__init__
++++++++++

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