import sys
if sys.version_info < (2,7) or not sys.version_info < (2,8):
    sys.exit('cnstools requires python 2.7')

import __init__ as cnstools
import _utils
import argparse

def main():
    program_dict = {}
    for module in cnstools.__all__:
        program_dict[module] = getattr(cnstools,module)

    main_parser = argparse.ArgumentParser(description="The cnstools package can be run in three ways. It can be run as an executable directly from the commandline, run as a python program with ``$python cnstools`` or it can be imported as a python module :class:`cnstools` and used with other python scripts.")
    prog_sub = main_parser.add_subparsers(title="Program",dest="program",help="Specifies which program in the cnstools package should be run.")

    program_run_dict =    {name:program_dict[name]._cl_entry for name in program_dict}
    program_parser_dict = {name:program_dict[name]._parser(prog_sub.add_parser,name) for name in program_dict}

    args = main_parser.parse_args(sys.argv[1:])

    #_utils.header_print("Running %s..."%args.program,h_type=3)

    program_to_run = args.program
    del args.program
    arg_dict = vars(args)
    arg_dict_to_pass = {key:arg_dict[key] for key in arg_dict if key!="program" and arg_dict[key]!=None}
    program_run_dict[program_to_run](**arg_dict_to_pass)

if __name__ == '__main__':
    main()
    sys.stderr.flush()
    sys.stdout.flush()