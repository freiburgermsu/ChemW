from pandas import DataFrame
from shutil import rmtree
from math import isclose
from glob import glob
#import requests, io, json
import chemw
import re, os


def test_inits():
    # import the class modules
    chem_mw = chemw.ChemMW()
    phreeq_db = chemw.PHREEQdb()
    print(os.getcwd())
    for TF in [chem_mw.verbose, chem_mw.final, chem_mw.end, phreeq_db.verbose]:
        assert type(TF) is bool
        
    rmtree(phreeq_db.output_path)
        

def test_accuracy():
    # calculate the MW for chemicals of known MW 
    test_chemicals = {
        'Na2.43_Cl_(OH)2_(OH)1.2_(OH)': 162.7,
        'Na2.43Cl(Ca(OH)2)1.2':180.2,
        'Na2.43Cl:2H2O': 127.3,
        'Na2.43Cl2.5:2H2O': 180.5,
        'CaCl2:(MgCl2)2:12H2O': 517.6,
        'Na2SO4:3K2SO4': 664.8,
        'K2SO4:CaSO4:H2O': 328.4,
        'Na.96Al.96Si2.04O6:H2O ': 219.2,
        'Ca1.019Na.136K.006Al2.18Si6.82O18:7.33H2O': 714.4
    }
    
    # calculate the MW for the dictionary of chemicals    
    chem_mw = chemw.ChemMW()
    for chemical in test_chemicals:
        chem_mw.mass(chemical)
        tolerance = chem_mw.mw*0.001 # 99.9% accuracy
        if not isclose(chem_mw.raw_mw, test_chemicals[chemical], rel_tol = tolerance):
            assert False
        else:
            assert True
            
    # affirm that iterated entities are zero
    for zero in [chem_mw.groups, chem_mw.layer, chem_mw.skip_characters]:
        assert zero == 0
        
    # test the MW of common names
    common_chemicals = {
    'water': 18.01528, 
    'acetone': 58.07914, 
    'toluene': 92.13842, 
    'glucose': 180.15588, 
    'sucrose': 342.29648, 
    'aspirin': 180.15742, 
    'hydrochloric acid': 36.46094,
    "alanine": 89.09318,
    "arginine": 174.20096,
    "asparagine": 132.11792,
    "aspartic acid": 133.10268,
    "cysteine": 121.15818,
    "glutamic acid": 147.12926,
    "glutamine": 146.1445,
    "glycine": 75.0666,
    "histidine": 155.15456,
    "isoleucine": 131.17292,
    "leucine": 131.17292,
    "lysine": 146.18756,
    "methionine": 149.21134,
    "phenylalanine": 165.18914,
    "proline": 115.13046,
    "serine": 105.09258,
    "threonine": 119.11916,
    "tryptophan": 204.22518,
    "tyrosine": 181.18854,
    "valine": 117.14634,
    }
    for chem in common_chemicals:
        assert chem_mw.mass(common_name = chem) == common_chemicals[chem]
            
            
def test_phreeq_db():
    # process the PHREEQ databases 
    phreeq_databases = [db for db in glob('databases/*.dat')]
    phreeq_db = chemw.PHREEQdb()
    for db in phreeq_databases:
        print('\n\n\n', re.search('([A-Za-z0-9_\.]+(?=\.dat))',db).group(), 'database\n', '='*len(db))
        phreeq_db.process(db)
        
    # verify the output folder and its contents 
    for db in phreeq_databases:
        json_name = re.search('([A-Za-z0-9_\.]+(?=\.dat))', db).group()+'.json'
        assert os.path.exists(os.path.join(phreeq_db.output_path, json_name))
        assert type(phreeq_db.db_name) is str
        assert type(phreeq_db.db) is DataFrame
   
    # delete the directory 
    rmtree(phreeq_db.output_path)
    
def test_proteins():
    protein = chemw.Proteins()
    
    # verify the sequence accuracy
    assert protein.mass('VPVIHTKPLFIRNFDQRCSCSRCFYLHSSTYIECTYISRFSKISLVSVTDFSLNGNVSTVFVPATRDSVPLHIIAPSSLIV') == 10537.67404
    
    # verify the FASTA accuracy
    protein.mass(fasta_path = os.path.join(os.path.dirname(__file__),'protein_sequence.fasta'))
    for line in protein.fasta_lines:
        if re.search('([0-9.]+)(?=_amu)', line):
            mass = float(re.search('([0-9.]+)(?=_amu)', line).group())
        else:
            seq = line.rstrip()
            assert protein.fasta_protein_masses[seq] == mass