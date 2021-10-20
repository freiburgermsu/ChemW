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