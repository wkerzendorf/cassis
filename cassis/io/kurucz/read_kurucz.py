import re
from astropy import units as u, constants as const
from cassis.pyparsing import Word, Literal, Group, Keyword, OneOrMore, Optional, Combine, alphas, nums, alphanums
from collections import OrderedDict
import pandas as pd


plusorminus = Literal( '+' ) | Literal( '-' )
decimal_pattern = Combine(Optional(plusorminus) + Word(nums, nums+",") + Optional("." + OneOrMore(Word(nums))))
decimal_pattern.setParseAction(lambda token: ''.join(token), lambda token: float(token[0]))

int_pattern = Word(nums).addParseAction(lambda token: int(token[0]))


atomic_number_pattern = Word(nums, min=1, max=3).setResultsName('atomic_number').setParseAction(lambda i: int(i[0]) )
ion_number_pattern = Word(nums, min=1, max=3).setResultsName('ion_number').setParseAction(lambda i: int(i[0]) )

species_pattern = atomic_number_pattern + Literal('.').suppress() + ion_number_pattern


header_pattern = species_pattern + Word(nums).setResultsName('lines_saved').setParseAction(lambda i: int(i[0])) + Literal('lines saved').suppress()
header_pattern += Word(nums).setResultsName('positive_lines_saved').setParseAction(lambda i: int(i[0])) + Literal('positive lines saved').suppress()
header_pattern += Word(nums).setResultsName('even_lines').setParseAction(lambda i: int(i[0])) + Literal('even').suppress()
header_pattern += Word(nums).setResultsName('odd_lines').setParseAction(lambda i: int(i[0])) + Literal('odd levels').suppress()
header_pattern += Word(nums).setResultsName('ionization_potential').setParseAction(lambda i: ((float(i[0])/u.cm) * const.c * const.h).to('eV')  )  + Keyword('ion pot cm-1 eve').suppress()

term_pattern = Group(OneOrMore(Word(alphanums, min=2))).setResultsName('term')
term_label_pattern = Word(alphanums, min=1, max=1).setResultsName('label') + term_pattern

#float_pattern = decimal_pattern

level_energy = decimal_pattern.copy()
level_energy.addParseAction(lambda token: (float(token[0])/u.cm * const.c * const.h).to('eV'))

level_data = species_pattern + (Literal('EVE') ^ Literal('ODD')).setResultsName('parity') 
level_data += int_pattern.setResultsName('index') 
level_data += decimal_pattern.setResultsName('energy') 
level_data += decimal_pattern.setResultsName('J')

def read_gf_gam_file(fname):
    """
    Reading a gfxxxx.gam file
    
    Parameters:
    -----------
    
    fname: str
        filename to read in 
        
    Returns:
    --------
    
    ~dict containing the header
    ~pandas.Dataframe object containing the tabular data
    
    
    """
    fhandle = open(fname)

    #reading the header
    for line in fhandle:
        header = header_pattern.parseString(line).asDict()
        break

    print 'header', header

    #reading the term label
    for line in fhandle:
        if line.strip().startswith('level'): break
#        print term_label_pattern.searchString(line)
        
    for line in fhandle:
        if line.strip().startswith('level'): continue
        if line.strip().startswith('ELEM'): break
    
    
    levels_data = OrderedDict([('atomic_number',[]) , ('ion_number',[]), ('parity',[]), ('index', []), ('energy', []), ('J', [])])
    
    
    for line in fhandle:
        level_data_dict = level_data.parseString(line).asDict()
        for key in level_data_dict:
            levels_data[key].append(level_data_dict[key])
    
    levels_data = pd.DataFrame(levels_data)
    levels_data['observed'] = levels_data['energy'] >= 0
    levels_data['energy'] = abs((levels_data['energy'].values/u.cm * const.c * const.h).to('eV')).value
    
    
    
    return header, levels_data
            
    
    


    