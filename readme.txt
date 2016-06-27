
Name: find_ncs.py
Requires: Python 2.7, progress_tracker.py

Usage: $.../find_ncs.py CHROMOSOME_NAME MAF_FILE GFF3_FILE OUTPUT_FILENAME
Usage (python 2.7):
    import find_ncs
    find_ncs.find_ncs(CHROMOSOME_NAME,MAF_FILE,GFF3_FILE,OUTPUT_FILENAME)

Description:
    Identifies conserved non-coding regions in the alignment data for a single chromosome (MAF_FILE) such that coding sequences identified in GFF3_FILE by the chromosome name CHROMOSOME_NAME are no longer included in any alignments. CNC ranges on the reference chromosome are output as inclusive ranges to the file OUTPUT_FILENAME. This file is tab delimited and in the format: start, stop, and index of parent sequence in the original .maf file. This file can be processed into a new .maf file using make_new_maf.py.

--------------------

Name: make_new_maf.py
Requires: Python 2.7, progress_tracker.py

Usage:
    $.../make_new_maf.py RANGE_LIST_FILE MAF_FILE OUTPUT_FILENAME
Usage (python 2.7):
    import make_new_maf
    make_new_maf.make_new_maf(RANGE_LIST_FILE,MAF_FILE,OUTPUT_FILENAME)

Description:
    Generate a new .maf file containing non-coding regions using RANGE_LIST_FILE (a find_ncs.py output file) and MAF_FILE (the .maf file originally used to generate that file using find_ncs.py).

--------------------

Name: progress_tracker.py
Requires: Python 2.7

Usage (python 2.7):
    import progress_tracker
    instance = progress_tracker.Progress_tracker(NAME,GOAL,[START_NOW])

Description:
    Contains a class for progress tracking and displaying it in terminal. Computes percent complete and gives a time estimate for loops. NAME will be displayed when stats are printed. GOAL is the number of steps it takes for a task to be complete. START_NOW begins timing of the first step immediately if True (defaults to False).

Methods:
    start()
    Description:
        Begins the timing of the first step

    step([NUMBER])
    Description:
        Tells the tracker that NUMBER steps have been completed (defaults to 1)

    display([SHOW_ESTIMATE])
    Description:
        Print name and progress (in percent (0-100) or "Finished") to the terminal and, if SHOW_ESTIMATE is True (as is default), print a time estimate as well.

    parseTime(seconds)
    Description:
        Mostly just an internal method, it could be repurposed. This is a simple recursive algorithm to parse seconds into a more human readable format.

--------------------
