Calculating the Molecular Weight from a Chemical Formula, Common Name, or Protein Sequence
-------------------------------------------------------------------------------------------------------------------------

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

The ``ChemMW`` object of ``ChemW`` parses a chemical formula string -- which can consist of any combination of elements and decimal stoichiometry -- and precisely calculates the MW of the chemical formula. The significant digits of the reported MW matches the lowest significant digits from the set of elements that constitute the chemical formula, where the used elemental masses are the most precise contemporary measurements of pure elements, per the `chemicals module <https://pypi.org/project/chemicals/>`_.

The ``Proteins`` and ``PHREEQdg`` objects are applications of the ``ChemMW`` object that expand the utility of this library. The ``Proteins`` object returns the mass of a protein by either parsing a string of a protein sequence, or by parsing a `FASTA-formatted file <https://en.wikipedia.org/wiki/FASTA_format>`_. This is applied in the `Codons module <https://pypi.org/project/codons/>`_ for genome-scale biology and bioengineering. The ``PHREEQdb`` object of ``ChemW`` parses a `PHREEQ database <https://www.usgs.gov/software/phreeqc-version-3>`_ via the ``ChemMW`` object and exports a JSON of mineral masses for all of the described minerals in the database. This is pivotally applied to calculate the masses of complex minerals in the PHREEQC databases of the `ROSSpy module <https://pypi.org/project/ROSSpy/>`_ for reverse osmosis research.

The ``ChemW`` library is offered with the `MIT License <https://opensource.org/licenses/MIT>`_\. Examples of the module are available in the examples directory of the `ChemW GitHub repository <https://github.com/freiburgermsu/ChemW>`_. Please submit errors or inaccuracies as `GitHub issues <https://github.com/freiburgermsu/ChemW/issues>`_ so that they may be resolved.


The full documentation is available on `ReadTheDocs <https://chemw.readthedocs.io/en/latest/index.html>`_.