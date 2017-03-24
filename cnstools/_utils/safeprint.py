import sys, threading

#printing with lock to make sure nothing overlaps/intertwines/interjects/etc
def safe_print(content,print_id=None):
    safe_print.print_lock.acquire()
    sys.stdout.flush()
    if  safe_print.print_id != print_id:
        safe_print.print_id = print_id
        if print_id!=None:
            sys.stdout.write("\n")
    sys.stdout.write(str(content)+"\n")
    sys.stdout.flush()
    safe_print.print_lock.release()
safe_print.print_id = None
safe_print.print_lock = threading.Lock()