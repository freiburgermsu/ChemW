from sigfig import round
from pandas import DataFrame
import chemw, re, os
from glob import glob
from shutil import rmtree


def test_inits():
    # import the class modules
    chem_mw = chemw.ChemMW()
    phreeq_db = chemw.PHREEQdb(r'..\examples\databases\pitzer.dat')
    
    # assert qualities of the modules
    assert type(phreeq_db.db_name) is str
    assert type(phreeq_db.db) is DataFrame
    for zero in [chem_mw.groups, chem_mw.layer, chem_mw.skip_characters]:
        assert zero == 0
    
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
    mw = chemw.ChemMW(verbose = False)
    for chemical in test_chemicals:
        mw.mass(chemical)
        if (round(mw.mineral_mass, 4) == test_chemicals[chemical]) or (test_chemicals[chemical]-amu_tolerance < mw.mineral_mass < test_chemicals[chemical]+amu_tolerance):
            assert True
        else:
            assert False
            
            
def test_phreeq_db():
    # process the PHREEQ databases 
    phreeq_databases = [db for db in glob(r'..\examples\databases\*.dat')]
    for db in phreeq_databases:
        print('\n\n\n', re.search('([A-Za-z0-9_\.]+(?=\.dat))',db).group(), 'database\n', '='*len(db))
        chemw.PHREEQdb(db)
        
    # verify the output folder and its contents
    export_path = os.path.join(os.getcwd(), f'PHREEQdb-0') 
    for db in phreeq_databases:
        assert os.path.exists(os.path.join(export_path, re.search('([A-Za-z0-9_\.]+(?=\.dat))', db).group()+'.json'))
   
    # delete the directory 
    rmtree(export_path)