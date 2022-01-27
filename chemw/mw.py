from pubchempy import get_compounds
from chemicals import periodic_table
from math import inf
from glob import glob
import requests, io
import pandas
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

        
class ChemMW():
    def __init__(self, verbose = False, printing = True):
        self.verbose = verbose
        self.printing = printing
        self.final = self.end = False

    def _parse_stoich(self,formula, ch_number):
        ch_no = ch_number  # ch_number is the first digit of the stoich float
        stoich = ''
        skips = 0
        if self.verbose:
            print('first ch', formula[ch_no])
        if formula[ch_no] == '.':
            stoich = '0.'
            ch_no += 1
            if re.search('[a-z]', formula[ch_no-2]):
                skips -= 1
        if (self.final or formula[ch_no] == ':') and not isnumber(formula[ch_no]):
            stoich = 1
        else:
            while re.search('[0-9\.]', formula[ch_no]): 
                if self.verbose:
                    print('later ch', formula[ch_no])
                stoich += formula[ch_no]
                if ch_no == len(formula)-1:
                    self.final = True
                    break
                ch_no += 1

        skips += ch_no - ch_number
        if stoich == '':
            stoich = 1
        if self.verbose:
            print('stoich',stoich)
        stoich = float(stoich)
        if re.search('(\.0$)', str(stoich)):
            stoich = int(stoich)
        return skips, stoich
    
    def _final(self,formula,ch_no):
        if ch_no+1 == len(formula)-1:
            self.final = True
        if ch_no >= len(formula)-1:
            self.end = True
            
    def _significant_digits(self, mass):
        mass_sigfigs = len(re.sub('\.', '', str(mass)))
        self.sigfigs = min(mass_sigfigs, self.sigfigs)

    def _group_parsing(self,formula, ch_number):
        self.group_masses = {
            1: 0,
            2: 0,
            3: 0,
            4: 0,
            5: 0,
            6: 0
        }
        self.layer += 1
        self.groups += 1
        ch_no2 = ch_number+1
        skip_characters = final_mass = 0
        if self.verbose:
            print('\n------Group parsing------')
        while self.layer > 0:
            ch = formula[ch_no2]
            if skip_characters > 0:
                if self.verbose:
                    print('\nch_no2', ch_no2)
                    print('skipping_group', ch)
                skip_characters -= 1
                ch_no2 += 1
                continue
            if self.verbose:
                print('\nch_no2', ch_no2)
                print('ch_2', ch)
                print(f'Layer {self.layer} mass: {self.group_masses[self.layer]}')
            if re.search('[\(\s]', ch):
                ch_no2 += 1
                self.layer += 1
                self.groups += 1
                self._final(formula,ch_no2)
                
                # parse the new sub mineral
                skips, mass = self._element_parsing(formula, ch_no2)
                self.group_masses[self.layer] += mass
                skip_characters += skips
                if self.verbose:
                    print('here1')
                    print('skips', skips)
                    print('skip_characters_group', skip_characters)
                ch_no2 += 1
            elif re.search('[):]', ch):
                ch_no2 += 1
                stoich = 1
                self._final(formula,ch_no2)
                
                if not self.end:
                    skips, stoich = self._parse_stoich(formula, ch_no2)
                    if self.verbose:
                        print('skips', skips)
                    skip_characters += skips 
                if not self.final:
                    self.group_masses[self.layer] *= stoich
                if self.verbose:
                    print('here2')
                    print('skip_characters_group', skip_characters)
#                 if re.search('\)', ch):
                self.layer -= 1
            else:
                # set_trace()
                skips, mass = self._element_parsing(formula, ch_no2)
                self.group_masses[self.layer] += mass
                skip_characters += skips
                if self.verbose:
                    print('here3')
                    print('skips', skips)
                    print('skip_characters_group', skip_characters)
                ch_no2 += 1
               
            if ch_no2 >= len(formula)-1:
                self.end = True
                self.layer -= 1
                stoich = 1
#                 print(f'--> ERROR: The final character of {formula} is {formula[-1]}') 
        
        # calculate the final mass of the nested group
        for group in self.group_masses:
            if self.verbose:
                print(f'component {group} mass', self.group_masses[group])
            final_mass += self.group_masses[group]
        if not self.end:
            if not re.search('[A-Z:\(]', formula[ch_no2+1]):
                skips, stoich = self._parse_stoich(formula, ch_no2)
            else:
                stoich = 1
        else:
            skip_characters += 1
                
        if self.verbose:
            print('stoich', stoich)
            print('final_mass', final_mass)
        if self.group_masses[2] > 0:
            final_mass *= stoich
        
        # calculate the total mass for nested groups 
        skip_characters = ch_no2 - ch_number + skip_characters-1
        if self.verbose:
            print('skip_characters_group', skip_characters)
            print('\n------End Group------\n')
        return skip_characters, final_mass

    def _element_parsing(self,formula, ch_number):
        stoich = skips = 0
        # set_trace()
        if re.search('[ +)]',formula[ch_number]):
            skip_characters = 0
            mass = 0
            if self.verbose:
                print('skip_characters_mineral_formula1', skip_characters)
            return skip_characters, mass
        elif re.search('[A-Z]',formula[ch_number]):
            self._final(formula,ch_number)
            element = formula[ch_number]
            if self.end:
                element = formula[ch_number]
                stoich = 1
                
                self._significant_digits(elemental_masses[element])
                mass = stoich * elemental_masses[element] 
                # track the elemental proportion
                if element not in self.element_masses:
                    self.element_masses[element] = mass
                else:
                    self.element_masses[element] += mass
                if self.verbose:
                    print('skip_characters_mineral_formula2', 0)
                return 0, mass

            elif re.search('[a-z]', formula[ch_number+1]):
                element = formula[ch_number] + formula[ch_number+1]          
                if ch_number+1 != len(formula)-1:
                    if re.search('[0-9]', formula[ch_number+2]):
                        skips, stoich = self._parse_stoich(formula, ch_number+2)
                    elif re.search('[A-Z\(:]', formula[ch_number+2]):
                        stoich = 1
                    elif formula[ch_number+2] == '.':
                        skips, stoich = self._parse_stoich(formula, ch_number+2)
                        skips += len('.')
                    elif formula[ch_number+2] == ' ':
                        stoich = 1
                        skips += len(' ')
                    else:
                        print('--> ERROR: The mineral formula {} may be unpredictable.'.format(formula))
                else:
                    stoich = 1

                self._significant_digits(elemental_masses[element])
                mass = stoich * elemental_masses[element] 
                
                # track the elemental proportion
                if element not in self.element_masses:
                    self.element_masses[element] = mass
                else:
                    self.element_masses[element] += mass

                skip_characters = skips + 1
                if self.verbose:
                    print('skip_characters_mineral_formula3', skip_characters)

                return skip_characters, mass

            elif re.search('[0-9]', formula[ch_number+1]):
                skips, stoich = self._parse_stoich(formula, ch_number+1)
                self._significant_digits(elemental_masses[element])
                mass = stoich * elemental_masses[element] 
                
                # track the elemental proportion
                if element not in self.element_masses:
                    self.element_masses[element] = mass
                else:
                    self.element_masses[element] += mass

                if self.verbose:
                    print('skip_characters_mineral_formula4', skips)
                return skips, mass

            elif formula[ch_number+1] == '.':
                skips, stoich = self._parse_stoich(formula, ch_number+1)
                self._significant_digits(elemental_masses[element])
                mass = stoich * elemental_masses[element] 
                
                # track the elemental proportion
                if element not in self.element_masses:
                    self.element_masses[element] = mass
                else:
                    self.element_masses[element] += mass

                skip_characters = skips
                if self.verbose:
                    print('skip_characters_mineral_formula5', skip_characters)
                return skip_characters, mass

            elif re.search('[A-Z():+ ]', formula[ch_number+1]):
                self._significant_digits(elemental_masses[element])
                mass = elemental_masses[element] 
                
                # track the elemental proportion
                if element not in self.element_masses:
                    self.element_masses[element] = mass
                else:
                    self.element_masses[element] += mass

                if self.verbose:
                    print('here4', stoich, elemental_masses[element])
                    print('skip_characters_mineral_formula6', 0)
                return 0, mass

        elif re.search(':',formula[ch_number]):
            skips = space = back_space = 0
            stoich = 1
            if re.search('[0-9]', formula[ch_number+1]):
                skips, stoich = self._parse_stoich(formula, ch_number+1)
                
                if self.verbose:
                    print('\npost-: value', formula[ch_number+skips+1])
                skip_characters, mass = self._group_parsing(formula, ch_number+skips)
                skip_characters += skips
                group_mass = mass * stoich
                if self.verbose:
                    print('skip_characters_mineral_formula9', skip_characters)
                    if not self.end:
                        print('post-skipping value', formula[ch_number+skip_characters])
                return skip_characters, group_mass
            
            if re.search('[ ]',formula[ch_number+1+skips]):
                space = 1
                back_space = 1
#             elif re.search('[A-Z]',formula[ch_number+skips+1]):
                
            else:
                if self.verbose:
                    print(f'--> ERROR: The {formula} formula is not predictable.')
                return 0, 0

        elif formula[ch_number-1] == '.':
            return 0, 0

        elif re.search('[0-9]', formula[ch_number]):
            self._final(formula,ch_number)            
            if not self.end:
                skips, stoich = self._parse_stoich(formula, ch_number)
            if self.verbose:
                print('skip_characters_mineral_formula', skips)
            return skips, stoich

    def _reset(self,):
        self.groups = self.layer = self.skip_characters = 0
        self.end = False
        self.final = False
        

    def mass(self, 
             formula: str = None,   # The molecular formula of the chemical whose mass will be calculated
             common_name: str = None  # The common name of the chemical, as they are recognized by PubChem
             ):  
        if common_name is not None:
            formula = pubchempy.get_compounds(common_name, 'name')[0].molecular_formula
        
        self.groups = self.layer = self.skip_characters = self.raw_mw = self.mw = 0 
        self.sigfigs = inf
        self.formula = formula
        formula = re.sub('[_]', '', formula)
        self.element_masses = {}
        if self.verbose:
            print('\n\n\n', formula, '\n', '='*2*len(formula))
        self.final = False
        for ch_number in range(len(formula)):
            ch = formula[ch_number]
#             if re.search('[\(:]', ch):
#                 self.skip_characters = 0
            if self.skip_characters > 0:
                self.skip_characters -= 1
                if ch_number == len(formula):
                    self.skip_characters = 0
                if self.verbose:
                    print('\nskip_characters', self.skip_characters)
                    print('skipping_mass', ch)
                continue
            if self.verbose:
                print('\ntotal_mass', self.raw_mw)
                print('ch_number', ch_number)
                print('ch', ch)

            if ch_number == len(formula)-1:
                self.final = True
                
            if formula[ch_number] == '(':
                self.skip_characters, mass = self._group_parsing(formula, ch_number)
                self.raw_mw += mass 
            else:
                self.skip_characters, mass = self._element_parsing(formula, ch_number)
                self.raw_mw += mass

        self.mw = round(self.raw_mw, self.sigfigs)
        if self.printing:
            print('{} --- MW (amu): {}'.format(formula, self.mw))
            
        # normalize the elemental proportions
        self.proportions = {}
        total_mass = sum([self.element_masses[element] for element in self.element_masses])
        for element in self.element_masses:
            self.proportions[element] = self.element_masses[element]/total_mass
            
        # reset the class object values
        self._reset()
        
        return self.mw
    
    
    
class Proteins():
    def __init__(self, verbose = False, printing = True):
        self.chem_mw = ChemMW(verbose = verbose, printing = printing)
        self.verbose = verbose        
        self.printing = printing
        
        # load amino acid masses, which were previously calculated through ChemMW to expedite computational time
        masses_path = os.path.join(os.path.dirname(__file__), 'amino_acids_masses.json')
        self.amino_acid_masses = json.load(open(masses_path))
        
    def _significant_digits(self, mass):
        mass_sigfigs = len(re.sub('\.', '', str(mass)))
        self.sigfigs = min(mass_sigfigs, self.sigfigs)
        
    def mass(self,
                protein_sequence: str = None, # the sequence of either one_letter or hyphenated three_letter amino acid sequences
                fasta_path: str = None, # the file to a local FASTA file
                fasta_link: str = None  # providing the link to a FASTA file as a string
                ):       
        def protein_mass(amino_acids):
            protein_mass = 0
            for amino_acid in amino_acids:
                if not re.search('[a-z]',amino_acid, flags = re.IGNORECASE):
                    if amino_acid != '*':
                        print(f'--> ERROR: An unexpected character {amino_acid} was encountered in the protein sequence {amino_acids}.')
                    continue
                mass = self.amino_acid_masses[amino_acid]
                self._significant_digits(mass)
                self.raw_protein_mass += mass
                protein_mass += mass
            return protein_mass
                
        if fasta_path is not None:
            with open(fasta_path) as input:
                self.fasta_lines = input.readlines()   
        elif fasta_link is not None:
            sequence = requests.get(fasta_link).content
            self.fasta_lines = io.StringIO(sequence.decode('utf-8')).readlines()
        else:
            three_letter_remainder = re.sub('(\-\w{3})', '', protein_sequence, flags = re.IGNORECASE)
            one_letter_remainder = re.sub('(\w)', '', protein_sequence, flags = re.IGNORECASE)
        
        self.raw_protein_mass = 0
        self.sigfigs = inf
        if fasta_path or fasta_link:
            self.fasta_protein_masses = {}
            for line in self.fasta_lines:
                if not re.search('>', line):
                    line = line.rstrip()
                    mass = protein_mass(line)
                    self.fasta_protein_masses[line] = mass
        elif three_letter_remainder == '' or three_letter_remainder == '*':
            amino_acids = protein_sequence.split('-')
            protein_mass(amino_acids)
        elif one_letter_remainder == '' or one_letter_remainder == '*':                
            protein_mass(protein_sequence)
        else:
            raise ImportError(f'The protein sequence {protein_sequence} has a remainder of {one_letter_remainder}, and does not follow the accepted conventions.')
            
        if fasta_path or fasta_link:
            for protein in self.fasta_protein_masses:
                self.fasta_protein_masses[protein] = round(self.fasta_protein_masses[protein], self.sigfigs)                
                if self.printing:
                    string = ' - '.join(['>Protein', f'{len(protein)}residues', f'{self.fasta_protein_masses[protein]}amu', f'\n{protein}'])
                    print(string)
            
            return self.fasta_protein_masses                
        else:
            self.protein_mass = round(self.raw_protein_mass, self.sigfigs)
            if self.printing:
                string = ' - '.join(['>Protein', f'{len(protein_sequence)}residues', f'{self.protein_mass}amu', f'\n{protein_sequence}'])
                print(string)
            
            return self.protein_mass


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
        if not os.path.exists(self.output_path):
            os.mkdir(self.output_path)
        

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

                    if re.search('(^\w+\s*\d*$)',self.db.at[index, 'content']):
                        reactants = self.db.at[index+1, 'content'].split(' = ')[0]
                        if all('#' not in entity for entity in reactants):
                            formula = reactants.split('+')[0].strip()
                            name = self.db.at[index, 'content']
                            name = re.sub('(\s+\d*)', '', name)
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
        self.db_name = re.search('([A-Za-z0-9_\.]+(?=\.dat))', db_path).group()
        self.db = pandas.read_table(db_path, sep='\n')
        self._database_parsing()

        # add elements to the JSON
        database_json = {'elements': {}, 'minerals': {}}
        for index, element in self.elements.iterrows():
            database_json['elements'][element['elements']] = {'charge_specie': element['alk'], 'gfw_formula':element['gfw_formula'], 'element_gfw':element['element_gfw']}

        # add minerals to the JSON
        for index, mineral in self.minerals.iterrows():
            phase = mineral['phases']                     
            if re.search('phases|phase', phase, flags = re.IGNORECASE):
                continue
                        
            # calculate the chemical masses for each mineral 
            formula = mineral['formula']
            formula = re.sub('Cyanide|Cyanate', 'CN', formula)   
            
            database_json['minerals'][phase] = {}
            database_json['minerals'][phase]['formula'] = formula
            database_json['minerals'][phase]['mass'] = self.chem_mw.mass(formula)

        # export the JSON files
        with open(os.path.join(self.output_path, f'{self.db_name}.json'), 'w') as output:
            json.dump(database_json, output, indent = 4)