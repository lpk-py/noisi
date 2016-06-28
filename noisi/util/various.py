from __future__ import print_function
import types
import re

def _get_val(prompt):
    """
    Separate call to raw_input to facilitate testing of input routine
    """
    return raw_input(prompt)

#ToDo: Return json dictionary instead of output string.
def _get_from_user(type,prompt,defval=None):
    """
    Prompt user to enter value that converts to a certain type until she
    enters a value valid for conversion to that type.
    :param type: The type of variable that should be returned (types object)
    :param prompt: Prompt that the user will be shown
    :param defval: Default value
    :return: A parameter read from command line converted to type.
    """
    istype = False
    while istype == False:
        val =  _get_val(prompt+' ['+str(defval)+']:  ')
        if val == '': val = defval
        try:
            val = type(val)
            istype = True
        except TypeError:
            print('\nMust enter input that converts to %s.\n' %type)            
            pass
        except ValueError:
            print('\n%s is no valid input, input must convert to %s.\n' %(val,type))            
            pass
    return(val)
    
    


def atoi(text):
    return int(text) if text.isdigit() else text

def natural_keys(text):
    '''
    Return integer number keys from lists that are strings.
    By: http://stackoverflow.com/questions/5967500/how-to-correctly-sort-a-string-with-a-number-inside
    '''
    return [ atoi(c) for c in re.split('(\d+)', text) ]
