from sigfig import round
from pandas import DataFrame
import chemw
import re, os
from glob import glob
from shutil import rmtree


def test_inits():
    # import the class modules
    chem_mw = chemw.ChemMW()
    phreeq_db = chemw.PHREEQdb()
    print(os.getcwd())
    for TF in [chem_mw.verbose, chem_mw.final, chem_mw.end, phreeq_db.verbose]:
        assert type(TF) is bool
        

def test_accuracy():
    # calculate the MW for chemicals of known MW 
    amu_tolerance = 1
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
        if (round(chem_mw.mw, 4) == test_chemicals[chemical]) or (test_chemicals[chemical]-amu_tolerance < chem_mw.mw < test_chemicals[chemical]+amu_tolerance):
            assert True
        else:
            assert False
            
    # affirm that iterated entities are zero
    for zero in [chem_mw.groups, chem_mw.layer, chem_mw.skip_characters]:
        assert zero == 0
            
            
def test_phreeq_db():
    # process the PHREEQ databases 
    phreeq_databases = [db for db in glob('databases/*.dat')]
    # output_path = os.path.join(os.path.dirname(__file__), 'PHREEQqb')
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