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
  },
  version = '0.1.1',
  license = 'MIT',
  description = "Calculates the Molecular Weight from a string of an arbitrary chemical formula.", 
  long_description = readme,
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/chemw',   
  keywords = ['chemistry', 'math', 'mass', 'weight', 'PHREEQC', 'molecular', 'mineral', 'formula', 'calculate'],
  install_requires = ['chemicals', 'pandas']
)