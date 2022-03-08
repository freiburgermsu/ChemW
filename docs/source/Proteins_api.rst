The ``Proteins`` package is an application of the ``ChemMW`` package. The ``Proteins`` object returns the mass of a protein, with its appropriate significant digits, by either parsing a string of a protein sequence, or by parsing a `FASTA-formatted file <https://en.wikipedia.org/wiki/FASTA_format>`_. This is applied in the `Codons module <https://pypi.org/project/codons/>`_ for genome-scale biology and bioengineering. 

Proteins
++++++++++++++++++

This class determines the MW of a protein sequence, after the initial parameters are specified: 

.. code-block:: python

  import chemw
  protein_mass = chemw.Proteins(verbose = False, printing = True)

- *verbose* & *printing* ``bool``: specifies whether troubleshooting information or MW results will be printed, respectively.

**Returns** *protein_mass* ``float``: The mass of the parameterized protein sequence to the appropriate significant digitsthat are alloted by those of the elemental masses.

++++++++++++++++
mass()
++++++++++++++++

The parameterized data is fitted to the Hill equation, with the following arguments and their default values:

.. code-block:: python

 protein.mass(protein_sequence = None,  fasta_path = None, fasta_link = None  # providing the link to a FASTA file as a string = None)

- *protein_sequence* ``str``: The protein sequence, with any combination of upper-case or lower-case letters, for which the MW is desired. The acceptable formats for the formula are quite broad, which are exemplified in the following formulae:

===================================================  ===================================================================================
 Example sequence                                                Format option
===================================================  ===================================================================================
 ``'LFCTHGLERVVZCLWHKRCCSTRLKSLLLRGCABC*'``            A single string of the one-letter amino acid codes. A trailing "*" is acceptable. 
``'gly-gln-his-ala-arg-asn-phe-pro-thr'``                A sequence of three-letter amino acid codes must be delimited with hyphens.
===================================================  ===================================================================================

- *fasta_path* & *fasta_link* ``str``: The path and URL link, respectively, to a FASTA file that contains the sequence, or multiple sequences, of the protein(s) for which the MW is desired. Each sequence must commence with a ``>`` as the first character of the description line.


++++++++++++++++++++++++++
Accessible content
++++++++++++++++++++++++++
The ``Proteins`` object retains numerous components that are accessible to the user: 

- *protein_mass* & *raw_protein_mass* ``float``: The protein mass that is adjusted and unadjusted for the appropriate number of significant digits.
- *fasta_protein_masses* ``dict``: A dictionary of each sequence from processing a FASTA file, where the value is the corresponding sequence's mass.
- *amino_acid_masses* ``dict``: A dictionary of all natural amino acids, and their masses to the appropriate number of significant digits.
- *fasta_lines* ``list``: The raw list of lines that constitute the loaded FASTA file, which can be used for post-processing.
- *sigfigs* ``float``: The number of sigfigs that are defined for each protein.
- *chem_mw* ``ChemMW``: An instance of the ``ChemMW`` object is loaded, which allows the user to access the ``ChemMW`` module through the ``PHREEQdb`` module.