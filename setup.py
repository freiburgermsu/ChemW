# -*- coding: utf-8 -*-
from setuptools import setup, find_packages

with open('README.rst') as file:
    readme = file.read()

setup(
  name = 'ChemW',      
  package_dir = {'mw':'chemw'},
  packages = find_packages(),
  version = '0.0.12',
  license = 'MIT',
  description = "Calculate the Molecular Weight from an arbitrary chemical formula as a string, and process PHREEQC databases into programmable JSON files.", 
  long_description = readme,
  author = 'Andrew Freiburger',               
  author_email = 'andrewfreiburger@gmail.com',
  url = 'https://github.com/freiburgermsu/chemw',   
  keywords = ['chemistry', 'math', 'mass', 'weight', 'PHREEQC', 'molecular', 'mineral'],
  install_requires = ['chemicals', 'pandas']
)