from chemicals import periodic_table
from glob import glob
from numpy import nan
from pdb import set_trace
import pandas
import json
import re


# define elemental masses
elemental_masses = {}
for element in periodic_table:
    elemental_masses[element.symbol] = element.MW
    
class phreeq_db():
    def __init__(self, db_path):
        self.db_name = re.search('[A-Za-z]+(?=\.)', db_path).group()
        self.db = pandas.read_table(db_path, sep='\n')
        self.chem_mw = ChemMW()
        
        # parse the database for elements and minerals
        self._database_parsing()
        self._database_json_creation()
        

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
        self.elements = pandas.DataFrame(elements_rows)
        self.elements.fillna(' ')
    #     elements.columns = elements.iloc[0]
        self.elements.drop([0], inplace = True)
        for column in self.elements:
            nan_entries = 0
            alphanumeric_entries = 0
            for entry in self.elements[column]:
    #             if entry in [None, ' ', nan]:
                if entry is not None:
                    if re.search('[a-z]|[0-9]', entry, re.IGNORECASE) and entry is not None:
        #             if entry is str or entry is float or entry is int:
                        alphanumeric_entries += 1
                    else:
                        nan_entries += 1
                else:
                    nan_entries += 1
            if nan_entries > alphanumeric_entries and len(self.elements.columns) > 5:
                print('deleted column: ', column)
                del self.elements[column]

        print(self.elements)
        self.elements.columns = ['elements', 'species', 'alk', 'gfw_formula', 'element_gfw']
    #     self.elements.rename(columns = {'SOLUTION_MASTER_SPECIES':'elements'}, inplace = True)
        print(self.elements)
    #         self.elements = self.elements.iloc[pandas.RangeIndex(len(self.elements)).drop([x for x in range(4)])]
        elements_list = list(self.elements['elements'])

        # define the minerals content for the database
        self.minerals = pandas.DataFrame(minerals_rows)
        self.minerals.columns = self.minerals.iloc[0]
        self.minerals = self.minerals.drop(0)
        print(self.minerals)
        mineral_list = list(self.minerals['phases'])
        formula_list = list(self.minerals['formula'])      
        
        
    def _database_json_creation():
        database_json = {'elements': {}, 'minerals': {}}

        # create the elements JSON
        for index, element in self.elements.iterrows():
            database_json['elements'][element['elements']] = {'charge_specie': element['alk'], 'gfw_formula':element['gfw_formula'], 'element_gfw':element['element_gfw']}

        # create the minerals JSON
        for index, mineral in self.minerals.iterrows():
            phase = mineral['phases']
            formula = mineral['formula']
            formula = re.sub('Cyanide|Cyanate', 'CN', formula)
            database_json['minerals'][phase] = {}
                        
            if re.search('PHASES', mineral['phases']):
                continue
            
            # calculate the chemical masses for each mineral 
            database_json['minerals'][phase]['formula'] = formula
            database_json['minerals'][phase]['mass'] = self.chem_mw.mass(formula)

        # export the JSON files
        with open(f'{self.db_name}.json', 'w') as output:
            json.dump(database_json, output, indent = 4)
        
        
class ChemMW():
    def __init__(self, database = False):
        self.database = database
        self.final = False
        self.groups = self.layer = self.skip_characters = 0

    def _parse_stoich(self,formula, ch_number):
        ch_no = ch_number  # ch_number is the first digit of the stoich float
        stoich = ''
        print('first ch', formula[ch_no])
        if formula[ch_no] == '.':
            stoich = '0.'
        if self.final:
            stoich = 1
        else:
            while re.search('[0-9\.]', formula[ch_no]): 
                print('later ch', formula[ch_no])
                stoich += formula[ch_no]
                if ch_no == len(formula)-1:
                    self.final = True
                    break
                ch_no += 1

        skips = ch_no - ch_number
        stoich = float(stoich)
        if re.search('(\.0$)', str(stoich)):
            stoich = int(stoich)
        return skips, stoich

    def _group_parsing(self,formula, ch_number):
        self.layer += 1
        self.groups += 1
        group_masses = {}
        ch_no2 = ch_number+1
        skip_characters = final_mass2 = group_masses[self.layer] = 0
        print('\nGroup parsing')
        while self.layer > 0:
            ch = formula[ch_no2]
            print('\nch_no2', ch_no2)
            print('ch_2', ch)
            if skip_characters > 0:
                print('skipping_group', skip_characters)
                skip_characters -= 1
                ch_no2 += 1
                continue
            if ch_no2+1 == len(formula)-1:
                self.final = True
                
            if re.search('[\(\s]', ch):
                ch_no2 += 1
                self.layer += 1
                self.groups += 1
                
                # parse the new sub mineral
                skips, mass = self._parse_mineral_formula(formula, ch_no2)
                group_masses[self.layer] = mass
                skip_characters += skips
                print('skips', skips)
                print('skip_characters_group', skip_characters)
            elif ch == ')':
                ch_no2 += 1
                stoich = 1
                if not self.final:
                    skips, stoich = self._parse_stoich(formula, ch_no2)
                    print('skips', skips)
                    skip_characters += skips 
                group_masses[self.layer] *= stoich
                print('skip_characters_group', skip_characters)
                self.layer -= 1
            else:
                # set_trace()
                skips, mass = self._parse_mineral_formula(formula, ch_no2)
                group_masses[self.layer] += mass
                print('skips', skips)
                skip_characters += skips
                print('skip_characters_group', skip_characters)
                ch_no2 += 1
               
    
#                 print(f'--> ERROR: The final character of {formula} is {formula[-1]}') 
        stoich = 1
        if not self.final:
            if not re.search('[A-Z:\(]', formula[ch_no2+1]):
                skips, stoich = self._parse_stoich(formula, ch_no2)

        # calculate the total mass for nested groups 
        final_mass = 0
        for group in group_masses:
            final_mass += group_masses[group]
        skip_characters = ch_no2 - ch_number + skip_characters-1
        print('skip_characters_group', skip_characters)
        
        print('\nEnd Group\n')
        return skip_characters, final_mass

    def _parse_mineral_formula(self,formula, ch_number):
        stoich = skips = 0
        # set_trace()
        if re.search('[ +)]',formula[ch_number]):
            skip_characters = 0
            mass = 0
            print('skip_characters_mineral_formula', skip_characters)
            return skip_characters, mass
        elif re.search('[A-Z]',formula[ch_number]):
            if self.final:
                element = formula[ch_number]
                stoich = 1
                mass = stoich * elemental_masses[element] 
                skip_characters = 0
                print('skip_characters_mineral_formula', skip_characters)
                return skip_characters, mass

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
                        skips = 1
                    else:
                        print('--> ERROR: The mineral formula {} is unpredictable.'.format(formula))
                else:
                    stoich = 1

                mass = stoich * elemental_masses[element] 
                skip_characters = skips + 1
                print('skip_characters_mineral_formula', skip_characters)

                return skip_characters, mass

            elif re.search('[0-9]', formula[ch_number+1]):
                skips, stoich = parse_stoich(formula, ch_number+1)

                element = formula[ch_number]
                mass = stoich * elemental_masses[element] 
                skip_characters = skips
                print('skip_characters_mineral_formula', skip_characters)
                return skip_characters, mass

            elif formula[ch_number+1] == '.':
                skips, stoich = parse_stoich(formula, ch_number+1)

                element = formula[ch_number]
                mass = stoich * elemental_masses[element] 
                skip_characters = skips+1
                print('skip_characters_mineral_formula', skip_characters)
                return skip_characters, mass

            elif re.search('[A-Z():+ ]', formula[ch_number+1]):
                element = formula[ch_number]
                stoich = 1
                mass = stoich * elemental_masses[element] 
                skip_characters = 0
                print('skip_characters_mineral_formula', skip_characters)
                return skip_characters, mass


        elif re.search(':',formula[ch_number]):
            skips = stoich = space = back_space = 0
            if re.search('[0-9]', formula[ch_number+1]):
                skips, stoich = parse_stoich(formula, ch_number+1)
            if re.search('[( ]',formula[ch_number+1+skips]):
                space = 1
                back_space = 1
            if formula[ch_number+1+skips+space:ch_number+4+skips+space] == 'H2O':
                skip_characters = len('H2O')+skips+back_space
                water_mass = elemental_masses['H'] * 2 + elemental_masses['O']
                mass = float(stoich) * water_mass
                print('skip_characters_mineral_formula', skip_characters)
                return skip_characters, mass
            elif formula[ch_number+1+skips+space:ch_number+4+skips+space] == 'H\+':
                skip_characters = len('H+')+skips
                proton_mass = elemental_masses['H']
                mass = float(stoich) * proton_mass
                print('skip_characters_mineral_formula', skip_characters)
                return skip_characters, mass
            elif re.search('[A-Z]',formula[ch_number+skips+1]):
                skip_characers, mass = group_parsing(formula, ch_number+skips, final)
                group_mass = mass * stoich
                return skip_characers, group_mass
                print('skip_characters_mineral_formula', skip_characters)
            else:
                print(f'--> ERROR: The {formula} formula is not predictable.')
                print('skip_characters_mineral_formula', skip_characters)
                return 0, 0

        elif formula[ch_number-1] == '.':
            return 0, 0

        elif re.search('[0-9]', formula[ch_number]):
            skips, stoich = self._parse_stoich(formula, ch_number)
            print('skip_characters_mineral_formula', skips)
            return skips, stoich


    def mass(self, formula):
        skip_characters = self.mineral_mass = 0 

        for ch_number in range(len(formula)):
            ch = formula[ch_number]
            print('\ntotal_mass', mineral_mass)
            print('ch_number', ch_number)
            print('ch', ch)
            if self.skip_characters > 0:
                self.skip_characters -= 1
                print('skip_characters', self.skip_characters)
                continue

            if ch_number == len(formula)-1:
                self.final = True
                
            if formula[ch_number] == '(':
                self.skip_characters, mass = self._group_parsing(formula, ch_number)
                self.mineral_mass += mass 
#                 if self.layer == 2:
#                     skip_characters += 1 
#                     if mineral == 'Berthierine_ISGS':
#                         skip_characters += 3
# #                     if mineral == 'Glauconite':
# #                         skip_characters += 1

#                 if self.layer == 1:
#                     if mineral in ['Boltwoodite', 'Corkite']:
#                         skip_characters -= 1
#                     if mineral == 'Glauconite':
#                         skip_characters -= 5
#                     if mineral in ['Saponite_SapCa', 'Vermiculite_SO']:
#                         skip_characters -= 9
#                     print('second')

#                 if self.layer == 0:
#                     if mineral in ['Brochantite', 'Borax', 'Antlerite', 'Corkite', 'Kasolite', 'Phosgenite', 'Tsumebite', 'Artinite', 'Jaffeite'] or (mineral in ['Burkeite', 'Dawsonite'] and re.search('sit', db)):
#                         skip_characters -= 1
#                     if mineral == 'Berthierine_ISGS':
#                         skip_characters -= 10
#                     if mineral == 'Glauconite':
#                         skip_characters -= 5
#                     if mineral in ['Saponite_SapCa', 'SmectiteMX80', 'Vermiculite_SO']:
#                         skip_characters -= 9
#                     print('first')
#                 elif mineral == 'Berthierine_ISGS':
#                     skip_characters -= 4

            else:
                self.skip_characters, mass = self._parse_mineral_formula(formula, ch_number)
                self.mineral_mass += mass

#         if mineral in ['Hydroxyapatite']:
#             minerals[mineral]['mass'] += elemental_masses['O']

        print('\n{} mass: {}'.format(formula, self.mineral_mass))

            
            
            
            
            
            
            
            
            
            
            
from chemicals import periodic_table
import pandas
import json
import re


# define elemental masses
elemental_masses = {}
for element in periodic_table:
    elemental_masses[element.symbol] = element.MW

class CpdMassPkg():
    def __init__(self):
        pass

    def mineral_mass(self, mineral_formula, mineral_name = None):
        self.formula = mineral_formula
        skip_characters = total_mass = mass = 0
        first = True
        double = triple = False
        for ch_number in range(len(self.formula)):
            print('total_mass', total_mass)
            print('ch_number', ch_number)
            print('character', self.formula[ch_number])
            if skip_characters > 0:
                print('skip_characters', skip_characters)
                skip_characters -= 1
                continue

            final = False
            if ch_number == len(self.formula)-1:
                print('final')
                final = True
            if self.formula[ch_number] == '(':
                skip_characters, mass = self.group_parsing(self.formula, ch_number, final)
                total_mass += mass 
                if triple:
                    skip_characters += 1 
                if double:
                    if mineral_name in ['Boltwoodite', 'Corkite']:
                        skip_characters -= 1
                        print('second')
                    triple = True
                if first:
                    if mineral_name in ['Brochantite', 'Borax', 'Antlerite', 'Corkite', 'Kasolite', 'Phosgenite', 'Tsumebite', 'Artinite']:
                        skip_characters -= 1
                    first = False
                    double = True
                    print('first')
            else:
                print('element', self.formula[ch_number])
                skip_characters, mass = self.parse_formula(self.formula, ch_number, final)
                total_mass += mass

        print(f'{self.formula} mass: {total_mass}')
        return total_mass

    def parse_stoich(self, formula, ch_number):
        ch_no = ch_number
        stoich = ''
        skip_minus = False
        if formula[ch_no] == '.':
            stoich = '0.'
            ch_no += 1
            skip_minus = True
        while re.search('[0-9\.]', formula[ch_no]):                    
            stoich += formula[ch_no]
            if ch_no == len(formula)-1:
                break
    #             if not re.search('[0-9\.]', formula[ch_no+1]):
    #                 break
            ch_no += 1

        skips = ch_no - ch_number 
        if skip_minus:
            skips -= 1
        stoich = float(stoich)
        if re.search('(\.0$)', str(stoich)):
            stoich = int(stoich)
        return skips, stoich

    def group_parsing(self, formula, ch_number, final):
        group_mass = group_mass2 = skip_characters = final_mass2 = total_skips = total_loops = skips2 = 0
        ch_no2 = ch_number + 1
        double = False
        group = True
        while formula[ch_no2] != ')' and ch_no2 != len(formula)-1:
            total_loops += 1
            print('ch_number', ch_no2)
            print('character', formula[ch_no2])
            print('skips_character_internal', skip_characters)
            if skip_characters > 0:
                skip_characters -= 1
                ch_no2 += 1
                total_skips += 1
                print('total_skips', total_skips)
                continue
            if re.search('[\(\s]', formula[ch_no2]):
                ch_no2 += 1
                double = True
                total_skips += 1
                print('total_skips1', total_skips)
                continue

            if double:
                skips, mass = self.parse_formula(formula, ch_no2, final, group = group)
                group_mass2 += mass
                skip_characters += skips
                ch_no2 += 1
                total_skips += 1
        #                     print(ch_no2)
                if formula[ch_no2] == ')':
                    print('group_mass2', group_mass2)
                    skips, stoich = self.parse_stoich(formula, ch_no2 + 1)
                    print('group_stoich2', stoich)
                    final_mass2 = group_mass2 * stoich
                    skip_characters += skips 
                    ch_no2 += 1
                    double = False
                total_skips += skip_characters
                print('total_skips2', total_skips)

            else:
                print(formula[ch_no2])
                skips, mass = self.parse_formula(formula, ch_no2, final, group = group)
                group_mass += mass
                skip_characters += skips
                ch_no2 += 1
                total_skips += skip_characters
                print('total_skips3', total_skips)
    #                     print(ch_no2)

        print('total_loops', total_loops)
        print('group_mass', group_mass)
        print('total_skips', total_skips)
        if ch_no2 == len(formula)-1 or re.search('[A-Z:\(]', formula[ch_no2+1]):
            stoich = 1
        else:
            skips2, stoich = self.parse_stoich(formula, ch_no2 + 1)
        print('group_stoich', stoich)
        print('skips2', skips2)
        final_mass = group_mass * stoich + final_mass2
    #     skip_characters = 2 + total_skips + skips2
        new_skips = skips2
        if skips2 == 0:
            new_skips = skips
        skip_characters = ch_no2 - ch_number + skips2 
        print('skip_characters_group', skip_characters)
    #     if double:
    #         skip_characters -= 1
        return skip_characters, final_mass

    def parse_formula(self, formula, ch_number, final = False, group = False):
        stoich = 0
        skips = 0
        if re.search('\s|\+',formula[ch_number]):
            skip_characters = 0
            mass = 0
            return skip_characters, mass
        elif re.search('[A-Z]',formula[ch_number]):
            if final or len(formula) == 1:
                element = formula[ch_number]
                print('\n', element)
                stoich = 1

                print('elemental_mass', elemental_masses[element])
                print('stoich', stoich)
                mass = stoich * elemental_masses[element] 
                skip_characters = 0
                print('mass', mass)
                return skip_characters, mass

            elif re.search('[a-z]', formula[ch_number+1]):
                element = formula[ch_number] + formula[ch_number+1]
                print('\n', element)            
                if ch_number+1 != len(formula)-1:
                    if re.search('[0-9]', formula[ch_number+2]):
                        skips, stoich = self.parse_stoich(formula, ch_number+2) # float(formula[ch_number+2])
    #                     if group:
    #                         skips = 0
                    elif re.search('[A-Z\(:]', formula[ch_number+2]):
                        stoich = 1
                    elif formula[ch_number+2] == '.':
                        skips, stoich = self.parse_stoich(formula, ch_number+2) # float(formula[ch_number+2])
                        skips += len('.')
                    elif formula[ch_number+2] == ' ':
                        stoich = 1
                        skips = 1
                    else:
                        print('--> ERROR: The mineral formula {} is unpredictable.'.format(formula))
                else:
                    stoich = 1

                print('elemental_mass1', elemental_masses[element])
                print('stoich1', stoich)
                mass = stoich * elemental_masses[element] 
                skip_characters = skips + 1
                print('mass1', mass)

                print('skips1', skip_characters)

                return skip_characters, mass

            elif re.search('[0-9]', formula[ch_number+1]):
                skips, stoich = self.parse_stoich(formula, ch_number+1)

                element = formula[ch_number]
                print('\n', element)
                print('elemental_mass2', elemental_masses[element])
                print('stoich2', stoich)
                mass = stoich * elemental_masses[element] 
                skip_characters = skips
                print('mass2', mass)
                return skip_characters, mass

            elif formula[ch_number+1] == '.':
                skips, stoich = parse_stoich(formula, ch_number+1)

                element = formula[ch_number]
                print('\n', element)
                print('elemental_mass2', elemental_masses[element])
                print('stoich2', stoich)
                mass = stoich * elemental_masses[element] 
                skip_characters = skips+1
                print('mass2', mass)
                return skip_characters, mass

            elif re.search('[A-Z\(\):\+\s]', formula[ch_number+1]):
                element = formula[ch_number]
                print('\n', element)
                stoich = 1

                print('elemental_mass3', elemental_masses[element])
                print('stoich3', stoich)
                mass = stoich * elemental_masses[element] 
                skip_characters = 0
                print('mass3', mass)
    #                 print('skips', skip_characters)
                return skip_characters, mass

        elif re.search(':',formula[ch_number]):
            print('element', formula[ch_number])
            skips = stoich = space = 0
            if re.search('[0-9]', formula[ch_number+1]):
                skips, stoich = self.parse_stoich(formula, ch_number+1)
            if ch_number+1+skips == ' ':
                space = 1
            if formula[ch_number+1+skips+space:ch_number+4+skips+space] == 'H2O':
                print('subbed portion', formula[ch_number+1+skips+space:ch_number+4+skips+space])
                skip_characters = len('H2O')+skips
                water_mass = elemental_masses['H'] * 2 + elemental_masses['O']
                mass = float(stoich) * water_mass
                return skip_characters, mass
            elif re.search('[A-Z]',formula[ch_number+skips+1]):
                skip_characers, mass = self.group_parsing(formula, ch_number+skips, final)
                group_mass = mass * stoich
                return skip_characers, group_mass
            else:
                print(f'--> ERROR: The {formula} formula is not predictable.')
                return 0, 0

        elif formula[ch_number-1] == '.':
            return 0, 0

        elif re.search('[0-9]', formula[ch_number]):
            skips, stoich = self.parse_stoich(formula, ch_number)
            return skips, stoich