Calculating the Molecular Weight from a Chemical Formula, Common Name, or Protein Sequence
-------------------------------------------------------------------------------------------------------------------------

|PyPI version| |Actions Status| |docs| |Downloads| |License|

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
   
.. |docs| image:: https://readthedocs.org/projects/chemw/badge/?version=latest
   :target: https://chemw.readthedocs.io/en/latest/?badge=latest
   :alt: Documentation Status


The ``ChemW`` library contains three packages, each with a distinct use-case, which are detailed in the following sections.

The ``ChemMW`` package parses any chemical formula that adheres to `chemical conventions <https://en.wikipedia.org/wiki/Chemical_formula>`_ -- which can consist of any combination of elements and decimal stoichiometry -- and calculates the MW of the chemical formula with the most accurate elemental masses that are available, per the `chemicals module <https://pypi.org/project/chemicals/>`_. The MW is precisely constrained to the significant digits of constituent elements in the chemical formula.

The ``Proteins`` and ``PHREEQdg`` packages are applications of the ``ChemMW`` package. The ``Proteins`` package returns the mass of a protein by either parsing a string of a protein sequence, or by parsing a `FASTA-formatted file <https://en.wikipedia.org/wiki/FASTA_format>`_. This is applied in the `Codons module <https://pypi.org/project/codons/>`_ for genome-scale biology and bioengineering. The ``PHREEQdb`` package parses a `PHREEQ database <https://www.usgs.gov/software/phreeqc-version-3>`_, calculates the MW of each mineral in the database, and exports a JSON of all mineral masses from the database. This is pivotally applied in the PHREEQC databases of the `ROSSpy module <https://pypi.org/project/ROSSpy/>`_ for reverse osmosis research.

The ``ChemW`` library is offered with the `MIT License <https://opensource.org/licenses/MIT>`_\. Examples of the module are available in the examples directory of the `ChemW GitHub repository <https://github.com/freiburgermsu/ChemW>`_. Please submit errors or inaccuracies as `GitHub issues <https://github.com/freiburgermsu/ChemW/issues>`_ so that they may be resolved.

Installation
----------------

The following command installs ``ChemW`` in a command prompt/terminal environment::
 
 pip install chemw


The full documentation is available on `ReadTheDocs <https://chemw.readthedocs.io/en/latest/index.html>`_.