import time
import sys
import threading

def safe_print(content,print_id=None):
    sys.stdout.flush()
    if not hasattr(safe_print, "print_id"): 
        safe_print.print_id = print_id
    elif safe_print.print_id != print_id:
        print "\n"
        safe_print.print_id = print_id
    print str(content)
    sys.stdout.flush()

def header_print(header):
    wall = lambda char: char*70
    safe_print("\n%s\n %s \n%s" % (wall("="),header,wall("~")))

class Progress_tracker():
    def __init__(self,name,size):
        self.name = name
        self.size = size
        self.estimate_on = True
        self.status_message = None
        self.steps = 0
        self.progress = 0
        self.start_time = None
        self.started = False
        self.is_done = False

    def start(self):
        if not self.started:
            self.start_time = time.time()
            self.started = True
        return self

    def done(self):
        self.steps = self.size
        self.ab = False
        self.display()
        self.is_done = True
        return self

    def step(self,size=1):
        if not self.started:
            self.start() 
        self.steps+=size
        return self

    def estimate(self,estimate_on):
        self.estimate_on = estimate_on
        return self

    def auto_display(self,rate=1):
        if rate>0:
            self.ad = True
            self.ad_rate = rate
            self.ad_next = time.time()
            self._auto_display()
        else:
            self.ad = False
        return self

    def _auto_display(self):
        if not self.is_done:
            if self.started:
                    self.display()
            if self.ad and self.steps<=self.size:
                self.ad_next = self.ad_next+self.ad_rate
                new_thread = threading.Timer(self.ad_next-time.time(),self._auto_display)
                new_thread.daemon = True
                new_thread.start()

    def display(self):

        curr_time = time.time()
        if not self.started: self.start_time = curr_time
        time_elapsed = curr_time-self.start_time

        self.progress = float(self.steps)/float(self.size)

        precentComplete = int((self.progress*100))
        tRemain = self._parseTime(int(time_elapsed/self.progress*(1-self.progress))) if self.progress!=0 else self._parseTime(float('inf'))

        write_string = "\033[F\x1b[2K%s:" % self.name
        write_string_size = len(self.name)+1
        write_string+= " "*(8*4-write_string_size)
        write_string_size = 8*4

        if self.steps<self.size:
            write_string += '\x1b[33m%s%s\x1b[0m' % (precentComplete, '%')
            write_string_size += len(str(precentComplete))+1
            if self.estimate_on:
                write_string+= " "*(8*5-write_string_size)
                write_string_size = 8*5
                rem = '(%s remaining)' % (tRemain)
                write_string+= rem
                write_string_size += len(rem)
            if self.status_message:
                write_string+= " "*(8*8-write_string_size)
                write_string_size = 8*8
                write_string+= '\x1b[34m[%s]\x1b[0m' % self.status_message
                write_string_size+=len(self.status_message)+2
        else:
            write_string += '\x1b[32mFinished\x1b[0m'

        safe_print(write_string,self)
        return self

    def _parseTime(self,seconds):
        if seconds < 60: return str(seconds//1)+"s"
        elif seconds == float('inf'): return "???"
        elif seconds < 3600: return str(seconds//60)+"m."+self._parseTime(seconds%60)
        elif seconds < 86400: return str(seconds//3600)+"h."+self._parseTime(seconds%3600)
        else: return str(seconds//86400)+"d."+self._parseTime(seconds%86400)

    def status(self, status_message=None):
        if not status_message: 
            status_message=None
        self.status_message = status_message
        self.display(estimate=self.estimate_on)
        return self

