from pubchempy import get_compounds
from chemicals import periodic_table
from warnings import warn
from pprint import pprint
from math import inf
from glob import glob
import requests, io
import pandas
import sigfig 
import json, re, os


def isnumber(num):
    try:
        float(num)
        return True
    except:
        try:
            int(num)
            return True
        except:
            return False

# define elemental masses
elemental_masses = {}
for element in periodic_table:
    elemental_masses[element.symbol] = element.MW
    
# allows case insensitive dictionary searches
class CaseInsensitiveDict(dict):        # sourced from https://stackoverflow.com/questions/2082152/case-insensitive-dictionary
    @classmethod
    def _k(cls, key):
        return key.lower() if isinstance(key, str) else key

    def __init__(self, *args, **kwargs):
        super(CaseInsensitiveDict, self).__init__(*args, **kwargs)
        CaseInsensitiveDict._convert_keys()
        
    def __getitem__(self, key):
        return super(CaseInsensitiveDict, self).__getitem__(__class__._k(key))
    
    def __setitem__(self, key, value):
        super(CaseInsensitiveDict, self).__setitem__(__class__._k(key), value)
        
    def __delitem__(self, key):
        return super(CaseInsensitiveDict, self).__delitem__(__class__._k(key))
    
    def __contains__(self, key):
        return super(CaseInsensitiveDict, self).__contains__(__class__._k(key))
    
    def has_key(self, key):
        return super(CaseInsensitiveDict, self).has_key(__class__._k(key))
    
    def pop(self, key, *args, **kwargs):
        return super(CaseInsensitiveDict, self).pop(__class__._k(key), *args, **kwargs)
    
    def get(self, key, *args, **kwargs):
        return super(CaseInsensitiveDict, self).get(__class__._k(key), *args, **kwargs)
    
    def setdefault(self, key, *args, **kwargs):
        return super(CaseInsensitiveDict, self).setdefault(__class__._k(key), *args, **kwargs)
    
    def update(self, E=None, **F):
        super(CaseInsensitiveDict, self).update(__class__(E))
        super(CaseInsensitiveDict, self).update(__class__(**F))
        
    def _convert_keys(self):
        for k in list(keys()):
            v = super(CaseInsensitiveDict, self).pop(k)
            CaseInsensitiveDict.__setitem__(k, v)

def _final(formula,ch_no):
    if ch_no+1 == len(formula)-1:
        final = True
    if ch_no >= len(formula)-1:
        end = True
    return final, end

def _significant_digits(mass, sigfigs=None):
    sigfigs = sigfigs or inf
    mass_sigfigs = len(re.sub(r'\.', '', str(mass)))
    sigfigs = min(mass_sigfigs, sigfigs)
    return sigfig.round(mass, sigfigs, warn=False)

def _element_check(element):
    # catch and ignore erroneous elements
    if element not in elemental_masses:
        warn(f'The character {element} is not an element.')
        return False
    return True

def _parse_stoich(formula, ch_number, verbose):
    ch_no = ch_number  # ch_number is the first digit of the stoich float
    stoich = ''
    skips = 0
    if verbose:
        print('first ch', formula[ch_no])
    if formula[ch_no] == '.':
        stoich = '0.'
        ch_no += 1
        if re.search('[a-z]', formula[ch_no-2]):
            skips -= 1
    final, end = _final(formula, ch_no)
    if (final or formula[ch_no] == ':') and not isnumber(formula[ch_no]):
        stoich = 1
    else:
        while re.search('[0-9.]', formula[ch_no]):
            if verbose:
                print('later ch', formula[ch_no])
            stoich += formula[ch_no]
            if ch_no == len(formula)-1:
                final = True
                break
            ch_no += 1

    skips += ch_no - ch_number
    if stoich == '':
        stoich = 1
    if verbose:
        print('stoich',stoich)
    stoich = float(stoich)
    if re.search(r'(\.0$)', str(stoich)):
        stoich = int(stoich)
    return skips, stoich

        
class ChemMW:

    def __init__(self):
        pass

    @staticmethod
    def _group_parsing(formula, ch_number, groups, layer, verbose):
        group_masses = {
            1: 0,
            2: 0,
            3: 0,
            4: 0,
            5: 0,
            6: 0
        }
        layer += 1
        groups += 1
        ch_no2 = ch_number+1
        skip_characters = final_mass = 0
        if verbose:
            print('\n------Group parsing------')
        while layer > 0:
            ch = formula[ch_no2]
            if skip_characters > 0:
                if verbose:
                    print('\nch_no2', ch_no2)
                    print('skipping_group', ch)
                skip_characters -= 1
                ch_no2 += 1
                continue
            if verbose:
                print('\nch_no2', ch_no2)
                print('ch_2', ch)
                print(f'Layer {layer} mass: {group_masses[layer]}')
            if re.search('[( ]', ch):
                ch_no2 += 1
                layer += 1
                groups += 1
                final, end = _final(formula,ch_no2)

                # parse the new sub mineral
                skips, mass = ChemMW._element_parsing(formula, ch_no2)
                group_masses[layer] += mass
                skip_characters += skips
                if verbose:
                    print('here1')
                    print('skips', skips)
                    print('skip_characters_group', skip_characters)
                ch_no2 += 1
            elif re.search('[):]', ch):
                ch_no2 += 1
                stoich = 1
                final, end = _final(formula,ch_no2)
                
                if not end:
                    skips, stoich = _parse_stoich(formula, ch_no2)
                    if verbose:
                        print('skips', skips)
                    skip_characters += skips 
                if not final:
                    group_masses[layer] *= stoich
                if verbose:
                    print('here2')
                    print('skip_characters_group', skip_characters)
#                 if re.search('\)', ch):
                layer -= 1
            else:
                # set_trace()
                skips, mass = ChemMW._element_parsing(formula, ch_no2)
                group_masses[layer] += mass
                skip_characters += skips
                if verbose:
                    print('here3')
                    print('skips', skips)
                    print('skip_characters_group', skip_characters)
                ch_no2 += 1
               
            if ch_no2 >= len(formula)-1:
                end = True
                layer -= 1
                stoich = 1
#                 print(f'--> ERROR: The final character of {formula} is {formula[-1]}') 
        
        # calculate the final mass of the nested group
        for group in group_masses:
            if verbose:
                print(f'component {group} mass', group_masses[group])
            final_mass += group_masses[group]
        if not end:
            if not re.search('[A-Z:(]', formula[ch_no2+1]):
                skips, stoich = _parse_stoich(formula, ch_no2)
            else:
                stoich = 1
        else:
            skip_characters += 1
                
        if verbose:
            print('stoich', stoich)
            print('final_mass', final_mass)
        if group_masses[2] > 0:
            final_mass *= stoich
        
        # calculate the total mass for nested groups 
        skip_characters = ch_no2 - ch_number + skip_characters-1
        if verbose:
            print('skip_characters_group', skip_characters)
            print('\n------End Group------\n')
        return skip_characters, final_mass

    @staticmethod
    def _element_parsing(formula, ch_number, element_masses, end, verbose):
        stoich = skips = 0
        # set_trace()
        if formula[ch_number] in ['*', "", None]:
            warn(f'The character {formula[ch_number]} is not defined.')
            return 0,0
        if re.search('[ +)]',formula[ch_number]):
            skip_characters = 0
            mass = 0
            if verbose:
                print('skip_characters_mineral_formula1', skip_characters)
            return skip_characters, mass
        elif re.search('[A-Z]',formula[ch_number]):
            final, end = _final(formula,ch_number)
            element = formula[ch_number]
            if end:
                element = formula[ch_number]
                stoich = 1
                
                if not _element_check(element):
                    return 0,0
                _significant_digits(elemental_masses[element])
                mass = stoich * elemental_masses[element] 
                # track the elemental proportion
                if element not in element_masses:
                    element_masses[element] = mass
                else:
                    element_masses[element] += mass
                if verbose:
                    print('skip_characters_mineral_formula2', 0)
                return 0, mass

            elif re.search('[a-z]', formula[ch_number+1]):
                element = formula[ch_number] + formula[ch_number+1]   
                if not _element_check(element):
                    return 0,0
                
                if ch_number+1 != len(formula)-1:
                    if re.search('[0-9]', formula[ch_number+2]):
                        skips, stoich = _parse_stoich(formula, ch_number+2)
                    elif re.search('[A-Z(:]', formula[ch_number+2]):
                        stoich = 1
                    elif formula[ch_number+2] == '.':
                        skips, stoich = _parse_stoich(formula, ch_number+2)
                        skips += len('.')
                    elif formula[ch_number+2] == ' ':
                        stoich = 1
                        skips += len(' ')
                    else:
                        warn(f'The mineral formula {formula} may be unpredictable.')
                else:
                    stoich = 1

                _significant_digits(elemental_masses[element])
                mass = stoich * elemental_masses[element] 
                
                # track the elemental proportion
                if element not in element_masses:
                    element_masses[element] = mass
                else:
                    element_masses[element] += mass

                skip_characters = skips + 1
                if verbose:
                    print('skip_characters_mineral_formula3', skip_characters)

                return skip_characters, mass

            elif re.search('[0-9]', formula[ch_number+1]):
                if not _element_check(element):
                    return 0,0
                skips, stoich = _parse_stoich(formula, ch_number+1)
                _significant_digits(elemental_masses[element])
                mass = stoich * elemental_masses[element] 
                
                # track the elemental proportion
                if element not in element_masses:
                    element_masses[element] = mass
                else:
                    element_masses[element] += mass

                if verbose:
                    print('skip_characters_mineral_formula4', skips)
                return skips, mass

            elif formula[ch_number+1] == '.':
                if not _element_check(element):
                    return 0,0
                skips, stoich = _parse_stoich(formula, ch_number+1)
                _significant_digits(elemental_masses[element])
                mass = stoich * elemental_masses[element] 
                
                # track the elemental proportion
                if element not in element_masses:
                    element_masses[element] = mass
                else:
                    element_masses[element] += mass

                skip_characters = skips
                if verbose:
                    print('skip_characters_mineral_formula5', skip_characters)
                return skip_characters, mass

            elif re.search('[A-Z():+ ]', formula[ch_number+1]):
                if not _element_check(element):
                    return 0,0
                _significant_digits(elemental_masses[element])
                mass = elemental_masses[element] 
                
                # track the elemental proportion
                if element not in element_masses:
                    element_masses[element] = mass
                else:
                    element_masses[element] += mass

                if verbose:
                    print('here4', stoich, elemental_masses[element])
                    print('skip_characters_mineral_formula6', 0)
                return 0, mass
            else:
                warn(f'ElementError: The formula {formula} character {formula[ch_number]} is not defined and will be skipped.')

        elif re.search(':',formula[ch_number]):
            skips = space = back_space = 0
            stoich = 1
            if re.search('[0-9]', formula[ch_number+1]):
                skips, stoich = _parse_stoich(formula, ch_number+1)
                
                if verbose:
                    print('\npost-: value', formula[ch_number+skips+1])
                skip_characters, mass = ChemMW._group_parsing(formula, ch_number+skips, groups)
                skip_characters += skips
                group_mass = mass * stoich
                if verbose:
                    print('skip_characters_mineral_formula9', skip_characters)
                    if not end:
                        print('post-skipping value', formula[ch_number+skip_characters])
                return skip_characters, group_mass
            
            if re.search('[ ]',formula[ch_number+1+skips]):
                space = 1
                back_space = 1
#             elif re.search('[A-Z]',formula[ch_number+skips+1]):
                
            else:
                if verbose:
                    warn(f'The {formula} formula is not predictable.')
                return 0, 0

        elif formula[ch_number-1] == '.':
            return 0, 0

        elif re.search('[0-9]', formula[ch_number]):
            _final(formula,ch_number)            
            if not end:
                skips, stoich = _parse_stoich(formula, ch_number)
            if verbose:
                print('skip_characters_mineral_formula', skips)
            return skips, stoich

    @staticmethod
    def mass(formula: str = None,      # The molecular formula of the chemical whose mass will be calculated
             common_name: str = None,  # The common name of the chemical, as they are recognized by PubChem
             final:bool = False,
             end:bool = False,
             verbose:bool = False,
             printing:bool = True,
             sigfigs:float = None
             ):
        sigfigs = sigfigs or inf
        if formula:
            if re.search('[a-z]',str(formula[0])):
                common_name = formula
        if formula in ['*', '']:
            return None
        if common_name:
            try:
                formula = get_compounds(common_name, 'name')[0].molecular_formula
            except:
                raise ValueError(f'The {common_name} common name is recognized by PubChem, and cannot be calculated through ChemW.')
        
        groups = layer = skip_characters = raw_mw = atoms = 0
        mw = ''
        formula = formula
        formula = re.sub('[_]', '', formula)
        element_masses = {}
        if verbose:
            print('\n\n\n', formula, '\n', '='*2*len(formula))
        final = False
        for ch_number in range(len(formula)):
            ch = formula[ch_number]
#             if re.search('[(:]', ch):
#                 skip_characters = 0
            if skip_characters > 0:
                skip_characters -= 1
                if ch_number == len(formula):
                    skip_characters = 0
                if verbose:
                    print('\nskip_characters', skip_characters)
                    print('skipping_mass', ch)
                continue
            if verbose:
                print('\ntotal_mass', raw_mw)
                print('ch_number', ch_number)
                print('ch', ch)

            if ch_number == len(formula)-1:
                final = True
                
            if formula[ch_number] == '(':
                skip_characters, mass = ChemMW._group_parsing(formula, ch_number, groups)
                raw_mw += mass 
            else:
                skip_characters, mass = ChemMW._element_parsing(formula, ch_number)
                raw_mw += mass

        mw = sigfig.round(str(raw_mw), sigfigs, warn = False)
        if printing:
            if common_name is not None:
                print('{}({}) --- MW (amu): {}'.format(common_name, formula, mw))
            else:
                print('{} --- MW (amu): {}'.format(formula, mw))
            
        # normalize the elemental proportions
        
        proportions = {}
        total_mass = sum([element_masses[element] for element in element_masses])
        for element in element_masses:
            proportions[element] = element_masses[element]/total_mass
            atoms += 1
            
        return mw
    
    
    
class Proteins():
    def __init__(self, verbose = False, printing = True):
        self.chem_mw = ChemMW(verbose = verbose, printing = printing)
        self.verbose = verbose
        self.printing = printing
        
        # load amino acid masses, which were previously calculated through ChemMW to expedite computational time
        with open(os.path.join(os.path.dirname(__file__), 'amino_acids_masses.json')) as aa_masses:
            self.amino_acid_masses = CaseInsensitiveDict(json.load(aa_masses))
        
    def mass(self,
                protein_sequence: str = None, # the sequence of either one_letter or hyphenated three_letter amino acid sequences
                fasta_path: str = None, # the file to a local FASTA file
                fasta_link: str = None  # providing the link to a FASTA file as a string
                ):       
        def calc_protein_mass(amino_acids):
            raw_protein_mass = 0
            for amino_acid in amino_acids:
                if not re.search('[A-Za-z]',amino_acid):
                    if amino_acid != '*':
                        print(f'--> ERROR: An unexpected character {amino_acid} was encountered in the protein sequence {amino_acids}.')
                    continue
                mass = self.amino_acid_masses.get(amino_acid)
                raw_protein_mass += mass
                sigs = _significant_digits(mass)
            return sigfig.round(str(raw_protein_mass), sigs, warn=False)
               
        self.fasta = []
        if fasta_path is not None:
            with open(fasta_path) as input:
                self.fasta_lines = input.readlines()
        elif fasta_link is not None:
            sequence = requests.get(fasta_link).content
            self.fasta_lines = io.StringIO(sequence.decode('utf-8')).readlines()
        else:
            remainder = re.sub(r'(\w)', '', protein_sequence, flags = re.IGNORECASE)
        
        self.raw_protein_mass = 0
        if fasta_path or fasta_link:
            self.fasta_protein_masses = {}
            protein = ''
            for line in self.fasta_lines:
                if '>' not in line:
                    line = line.rstrip()
                    protein += line
                elif '>' in line:
                    self.fasta_protein_masses[protein] = calc_protein_mass(protein)
                    protein = ''
                    
            self.fasta_protein_masses[protein] = calc_protein_mass(protein)
            if '' in self.fasta_protein_masses:
                self.fasta_protein_masses.pop('')
        elif remainder == '' or remainder == '*':    
            protein_mass = calc_protein_mass(protein_sequence)
        elif re.search(r'(-)+(\*)?', remainder):
            amino_acids = protein_sequence.split('-')
            protein_mass = calc_protein_mass(amino_acids)
            if self.printing:
                string = ' - '.join(['>Protein', f'{len(amino_acids)}_residues', f'{protein_mass}_amu']) + f'\n{protein_sequence}'
                print(string)
                self.fasta.append(string)
            self.fasta = '\n'.join(self.fasta)
            return protein_mass
        else:
            raise ImportError(f'The protein sequence {protein_sequence} has a remainder of '
                              f'{one_letter_remainder}, and does not follow the accepted conventions.')
            
        if fasta_path or fasta_link:
            if self.printing:
                for protein in self.fasta_protein_masses:
                    string = ' - '.join(['>Protein', f'{len(protein)}_residues',
                                         f'{self.fasta_protein_masses[protein]}_amu']) + f'\n{protein}'
                    print(string)
                    self.fasta.append(string)
            self.fasta = '\n'.join(self.fasta)
            return self.fasta_protein_masses
        else:
            if self.printing:
                string = ' - '.join(['>Protein', f'{len(protein_sequence)}_residues', f'{protein_mass}_amu']) + f'\n{protein_sequence}'
                print(string)
                self.fasta.append(string)
            self.fasta = '\n'.join(self.fasta)
            return protein_mass
        
class PHREEQdb():
    def __init__(self, output_path = None, verbose = False, printing = False):
        self.chem_mw = ChemMW(verbose = verbose, printing = printing)
        self.verbose = verbose
        
        # define the output path
        if output_path is None:
            count = 0
            self.output_path = os.path.join(os.getcwd(), f'PHREEQdb-{count}')
            while os.path.exists(self.output_path):
                count += 1
                self.output_path = os.path.join(os.getcwd(), f'PHREEQdb-{count}')
        else:
            self.output_path = output_path
        if not os.path.exists(output_path):
            os.mkdir(output_path)
        

    def _database_parsing(self,):
        start_master = False
        if self.db.columns == ['SOLUTION_MASTER_SPECIES']:
            start_master = True
        self.db.columns = ['content']

        elements_rows = []
        minerals_rows = []
        elemental_parsing = False
        mineral_parsing = False
        for index, row in self.db.iterrows():
            if (re.search('SOLUTION_MASTER_SPECIES', row['content']) or start_master) and not elemental_parsing :
                while not re.search('SOLUTION_SPECIES', self.db.at[index, 'content']):
                    split_row = self.db.at[index, 'content'].split()
                    if all(not re.search('^#', entity) for entity in split_row):
                        elements_rows.append(split_row)
                    index+= 1
                elemental_parsing = True

            if re.search('PHASES', row['content']) and not mineral_parsing:
                loop = False
                while not re.search('PITZER|EXCHANGE_MASTER_SPECIES|SURFACE_MASTER_SPECIES', self.db.at[index, 'content']):
                    if not loop:
                        minerals_rows.append(['phases', 'formula'])
                        loop = True
                        index += 1
                        continue

                    if re.search(r'(^\w+\s*\d*$)',self.db.at[index, 'content']):
                        reactants = self.db.at[index+1, 'content'].split(' = ')[0]
                        if all('#' not in entity for entity in reactants):
                            formula = reactants.split('+')[0].strip()
                            name = self.db.at[index, 'content']
                            name = re.sub(r'(\s+\d*)', '', name)
                            minerals_rows.append([name, formula])
                    index+= 1

                    if index == len(self.db):
                        break
                mineral_parsing = True

        # define the elements content for the database
        self.minerals = pandas.DataFrame(minerals_rows)
        self.minerals.columns = self.minerals.iloc[0]
        self.minerals = self.minerals.drop(0)
        
        self.elements = pandas.DataFrame(elements_rows)
        self.elements.fillna(' ')
        self.elements.drop([0], inplace = True)
        for column in self.elements:
            nan_entries = 0
            alphanumeric_entries = 0
            for entry in self.elements[column]:
                if entry is not None:
                    if re.search('[a-z]|[0-9]', entry, re.IGNORECASE):
                        alphanumeric_entries += 1
                    else:
                        nan_entries += 1
                else:
                    nan_entries += 1
            if nan_entries > alphanumeric_entries and len(self.elements.columns) > 5:
                print('deleted column: ', column)
                del self.elements[column]

        self.elements.columns = ['elements', 'species', 'alk', 'gfw_formula', 'element_gfw']
        
        if self.verbose:
            print(self.elements)
            print(self.minerals)
        
        
    def process(self,db_path):
        # load the database
        db_name = re.search(r'([A-Za-z0-9_.]+(?=\.dat))', db_path).group()
        with open(db_path, 'r', encoding="utf8") as db:
            db = db.readlines()
            self.db = pandas.DataFrame(db)
        self._database_parsing()

        # add elements to the JSON
        database_json = {'elements': {}, 'minerals': {}}
        for index, element in self.elements.iterrows():
            database_json['elements'][element['elements']] = {'charge_specie': element['alk'], 'gfw_formula':element['gfw_formula'], 'element_gfw':element['element_gfw']}

        # add minerals to the JSON
        minerals = set()
        for index, mineral in self.minerals.iterrows():
            phase = mineral['phases']                     
            if re.search('phases|phase', phase, flags = re.IGNORECASE):
                continue
                        
            # calculate the chemical masses for each mineral 
            formula = mineral['formula']
            formula = re.sub('Cyanide|Cyanate', 'CN', formula)   
            minerals.add(formula)
            
            database_json['minerals'][phase] = {}
            database_json['minerals'][phase]['formula'] = formula
            database_json['minerals'][phase]['mass'] = self.chem_mw.mass(formula)

        # export the JSON files
        print(db_name, f': {len(minerals)} minerals')
        with open(os.path.join(self.output_path, f'{db_name}.json'), 'w') as output:
            json.dump(database_json, output, indent=4)
            
        return minerals