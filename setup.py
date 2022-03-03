# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.rst') as file:
    readme = file.read()

setup(
  name = 'ChemW',      
  package_dir = {'mw':'chemw'},
  packages = find_packages(),
  package_data = {
	'test':['databases/*', 'protein_sequence.fasta'],
    'chemw':['amino_acids_masses.json'],
  },
  version = '0.2.3',
  license = 'MIT',
  description = "Calculates the Molecular Weight, to the appropriate significant digits, from a string of an arbitrary chemical formula, a protein sequence of one- or three-letter codes, or chemical common names that are recognized by PubChem.", 
  long_description = readme,
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/chemw',   
  keywords = ['chemistry', 'math', 'mass', 'weight', 'PHREEQC', 'molecular', 'mineral', 'formula', 'calculate'],
  install_requires = ['chemicals', 'pandas', 'pubchempy', 'requests']
)