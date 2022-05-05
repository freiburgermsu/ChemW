The ``ChemMW`` package of ``ChemW`` parses a chemical formula string -- which can consist of any combination of elements and decimal stoichiometry -- and precisely calculates the MW of the chemical formula. The significant digits of the reported MW matches the lowest significant digits from the set of elements that constitute the chemical formula, where the elemental masses are the most precise contemporary measurements of pure elements, per the `chemicals module <https://pypi.org/project/chemicals/>`_. This approach offers 1-2 more significant digits than `PubChemPy <https://pypi.org/project/PubChemPy/>`_).

ChemMW
++++++++++++++++++

This class determines the MW of a formula or common-name chemical, after the initial parameters are specified: 

.. code-block:: python

  import chemw
  chem_mw = chemw.ChemMW(verbose = False, printing = True)

- *verbose* & *printing* ``bool``: specifies whether troubleshooting information or MW results will be printed, respectively.

++++++++++++++++
mass()
++++++++++++++++

The parameterized data is fitted to the Hill equation, with the following arguments and their default values:

.. code-block:: python

 mw = chem_mw.mass(formula = None, common_name = None)

- *formula* ``str``: parameterizes the chemical formula for which the MW is desired. The acceptable formats for the formula are quite broad, which are exemplified in the following table:

===================================================  =========================================================================================================
 Example chemical                                      Format option
===================================================  =========================================================================================================
 ``'C6H12O6'``                                           Any organic compounds can be easily processed.
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
===================================================  =========================================================================================================
                                                            
- *common_name* ``str``: parameterizes the common name of the chemical for which the MW is desired, as it is recognized by `Pubchem <https://pubchem.ncbi.nlm.nih.gov>`_: e.g. ``'water'``, ``'acetone'``, ``'toluene'``, ``'glucose'``, ``'aspirin'``, ``'hydrochloric acid'``, ``"alanine"``, ``"glutamine"``, ``"phenylalanine"``, ``"tryptophan"``, et cetera.

**Returns** *mw* ``float``: The mass of the parameterized chemical to the appropriate significant digitsthat are alloted by those of the elemental masses.


++++++++++++++++++++++++++
Accessible content
++++++++++++++++++++++++++
The ``ChemMW`` object retains numerous components that are accessible to the user: 

- *mw* ``str``: The MW of the parameterized chemical formula with the proper significant digits, which should be used for printing and storing a calculated mass.
- *raw_mw* ``float``: The MW of the parameterized chemical formula without consideration of the significant digits, which should be used for calculations.
- *proportions* ``dict``: The ratio of elements in the chemical formula. This loses accuracy with the grouped elements, and is being improved.
- *formula* ``str``: The original chemical formula as a string.
- *atoms* ``int``: The qantity of atoms that are defined in the chemical formula.
- *groups* ``int``: A numerical counter for the quantity of chemical groups that are 
- *group_masses* ``dict``: A dictionary for the masses of each nesting level in a molecule.