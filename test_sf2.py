#!/usr/bin/env python 

'''Tests the installation and searches for dependencies in the path'''

#Imports
import imp
import os

#Functions
def check_basic():
    '''Assures the user if this script is run'''
    print '\n\033[0;1mStructureFold2 is in the path with proper permissions.\033[0m\n'

def py_module_status(amodule):
    '''Checks for the presence of modules'''
    info = {True:'The {} Python module is installed properly.',
            False:'The {} Python module is either not installed or misconfigured.'}
    try:
        imp.find_module(amodule)
        present = True
    except ImportError:
        present = False
    print info[present].format(amodule)

def check_py_list_modules(alist):
    '''runs py_module_status on each element in the list'''
    print '\033[0;1mTesting for Python dependencies...\033[0m'
    for item in alist:
        py_module_status(item)
    print ''

def check_dependencies(prog_dictionary):
    '''runs program_status on each element in the list'''
    print '\033[0;1mTesting for program dependencies...\033[0m'
    for program, error_message in prog_dictionary.items():
        print test(program,error_message)
    print ''

def is_exe(fpath):
    '''tests'''
    return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

def test(program,message):
    '''tests for the programs in the path'''
    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
            return '{} is located at {}'.format(program,exe_file)
    return 'Cannot find a valid install of {}, see {} or check permissions'.format(program,message)

def check_optional(options_dictionary):
    '''runs program_status on each element in the list'''
    print '\033[0;1mTesting for optional programs...\033[0m'
    for program, error_message in options_dictionary.items():
        print test_option(program,error_message)
    print ''

def test_option(program,message):
    '''tests for the programs in the path'''
    for path in os.environ["PATH"].split(os.pathsep):
        exe_file = os.path.join(path, program)
        if is_exe(exe_file):
            return '{} is located at {}'.format(program,exe_file)
    return 'Cannot find a valid install of {}, see {} or check permissions'.format(program,message)


def main():
    check_basic()
    modules = ['numpy','Bio']
    check_py_list_modules(modules)
    programs = {'bowtie2':'http://bowtie-bio.sourceforge.net/bowtie2/index.shtml',
                'cutadapt':'http://cutadapt.readthedocs.io/en/stable/guide.html',
                'samtools':'http://samtools.sourceforge.net/'}
    check_dependencies(programs)
    options = {'Fold':'https://rna.urmc.rochester.edu/RNAstructure.html'}
    check_optional(options)

if __name__ == '__main__':
    main()
