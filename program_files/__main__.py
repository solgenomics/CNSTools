import sys

import cnstools
program_run_dict = {}
for sub in cnstools.__all__:
    exec "program_run_dict['"+sub+"'] = cnstools."+sub+".run"

def main():
  program_run_dict[sys.argv[1].strip()](sys.argv[1:])
  
if __name__ == '__main__':
  main()