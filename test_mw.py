from sigfig import round
import mw

test_chemicals = {
    'Na2.43_Cl_(OH)2_(OH)1.2_(OH)': 162.7,
    'Na2.43Cl(Ca(OH)2)1.2':180.2,
    'Na2.43Cl:2H2O': 127.3,
    'Na2.43Cl2.5:2H2O': 180.5,
}


def test_mass():
    mw = ChemMW(verbose = False)
    for chemical in test_chemicals:
        print('\n\n\n', chemical, '\n', '='*2*len(chemical))
        mw.mass(chemical)
        if round(mw.mineral_mass, 4) != test_chemicals[chemical]:
            mw = ChemMW(verbose = True)
            mw.mass(chemical)
        else:
            print('pass\n\n')