#David Lyon    
                                                                                                                                                                                 
import time
import sys

last_displayed_progress_tracker = [None]

class Progress_tracker():
    """docstring for Progress"""
    def __init__(self,name,size,startIt=True):
        global last_displayed_progress_tracker
        self.last_displayed_progress_tracker = last_displayed_progress_tracker
        self.name = name
        self.size = size
        self.estimate_on = True
        self.prevStep = None
        self.status_message = None
        self.timeElapsed = 0
        self.steps = 0
        self.prog = 0.000001
        self.prevDisplay = 0
        self.autoDisplay = False
        self.displayRate = 0
        if startIt: self.start()
    def start(self):
        self.prevStep = time.time()
        return self
    def step(self,size=1):
        if not self.prevStep: self.start()
        elif self.steps<self.size:
            thisStep = time.time()
            self.timeElapsed+= thisStep-self.prevStep
            self.prevStep = thisStep
            self.steps += size
            if self.steps>self.size: self.steps=self.size
            self.prog = (self.steps/float(self.size))
            if self.autoDisplay and thisStep-self.prevDisplay>self.displayRate:
                self.display(finish=False)
        return self
    def display(self, estimate=True, rate=None, finish=True):
        if self.estimate_on!=estimate: self.estimate_on=estimate
        if rate!=None:
            if rate>0:
                self.autoDisplay = True
                self.displayRate = rate
            elif rate:
                self.autoDisplay = False
                self.displayRate = rate
        precentComplete = int((self.prog*100))
        tRemain = self.parseTime(int(self.timeElapsed/self.prog*(1-self.prog)))
        if self.last_displayed_progress_tracker[0]!=self: 
            self.last_displayed_progress_tracker[0] = self
            sys.stdout.write('\n')
            sys.stdout.flush()

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
            if finish: self.last_displayed_progress_tracker[0]=None
        print write_string
        sys.stdout.flush()
        self.prevDisplay = time.time()
        return self
    def parseTime(self,seconds):
        if seconds < 60: return str(seconds//1)+"s"
        elif seconds < 3600: return str(seconds//60)+"m."+self.parseTime(seconds%60)
        elif seconds < 86400: return str(seconds//3600)+"h."+self.parseTime(seconds%3600)
        else: return str(seconds//86400)+"d."+self.parseTime(seconds%86400)
    def status(self, status_message=None):
        if not status_message: status_message=None
        self.status_message = status_message
        self.display(estimate=self.estimate_on)
        return self




