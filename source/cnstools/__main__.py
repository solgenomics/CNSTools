import __init__ as cnstools
import sys
import _utils
import argparse

def main():

    program_dict = {}
    for module in cnstools.__all__:
        program_dict[module] = getattr(cnstools,module)
    program_run_dict = {name:program_dict[name].run for name in program_dict}
    program_parser_dict = {name:program_dict[name].parser for name in program_dict if hasattr(program_dict[name],"parser")}

    main_parser = argparse.ArgumentParser()
    main_parser.add_argument("program", help="program to run",choices=cnstools.__all__)

    main_args,program_arg_list = main_parser.parse_known_args(sys.argv[1:])

    _utils.header_print("Running %s..."%main_args.program)

    if main_args.program in program_parser_dict:
        program_args = program_parser_dict[main_args.program].parse_args(program_arg_list)
        program_run_dict[main_args.program](**vars(program_args))
    else:
        program_run_dict[main_args.program](*program_arg_list)

if __name__ == '__main__':
    main()
    sys.stderr.flush()
    sys.stdout.flush()