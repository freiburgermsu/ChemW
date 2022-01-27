# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.rst') as file:
    readme = file.read()

setup(
  name = 'ChemW',      
  package_dir = {'mw':'chemw'},
  packages = find_packages(),
  package_data = {
	'test':['databases/*'],
    'chemw':['amino_acids_masses.json'],
  },
  version = '0.1.5',
  license = 'MIT',
  description = "Calculates the Molecular Weight, to the appropriate significant digits, from a string of an arbitrary chemical formula, protein sequence, or common chemical name.", 
  long_description = readme,
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/chemw',   
  keywords = ['chemistry', 'math', 'mass', 'weight', 'PHREEQC', 'molecular', 'mineral', 'formula', 'calculate'],
  install_requires = ['chemicals', 'pandas', 'pubchempy', 'requests']
)