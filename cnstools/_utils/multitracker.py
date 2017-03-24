import threading
import time
from multiprinter import MultiPrinter

class _Tracker(object):
    """docstring for MultiTracker"""
    styles = ["percent","fraction","noProg"]
    spinner = ["...","..",".",".."]
    red, green, yellow, reset = "\033[0;31m","\033[0;32m","\033[0;33m","\033[0;m"#"","","",""
    def __init__(self,name,size,parent=None,estimate=True,style="percent",lock=threading.RLock()):
        with lock:
            self.name = name
            self.size = size if size!=0 else 1
            self.steps = 0
            self.parent = parent
            self._estimate_default = estimate
            self._estimate = False
            self.running = False
            self.start_time = None
            self.pause_adjustment = 0 #stores the previously elapsed time when the tracker is paused.
            self._done = False
            self._auto_display = False
            self.spinner_step = 0 if self.parent == None else self.parent.spinner_step
            self.children = []
            self.style = style.strip().split()
            self.lock = lock

    def toggle(self,val=None):
        with self.lock:
            if val == None: val = (not self.running)

            if val:
                self.start_time = time.time() - self.pause_adjustment
                self.pause_adjustment = 0
                self.running = True
            else:
                self.pause_adjustment = time.time() - self.start_time
                self.running = False

            return self

    def done(self):
        with self.lock:
            self._done = True
            self.steps = self.size
            self.children = []
            return self

    def step(self,size = 1):
        with self.lock:
            if not self.running: 
                self.toggle() 
                self._estimate = self._estimate_default
            self.steps += size
            if self.steps>=self.size: self.done()
            return self

    def estimate(self,val):
        with self.lock:
            self._estimate = val
            return self

    def __str__(self):
        with self.lock:
            res = self.name
            ex = 0
            self.spinner_step=(self.spinner_step+1)%len(self.spinner)
            if self._done:
                res+=":"
                res+=" "*max(1,30-len(res)+ex)
                res+=self.green+"Done!"+self.reset
                ex += len(self.green+self.reset)
            elif "percent" in self.style:
                res+=":"
                res+=" "*max(1,30-len(res)+ex)
                precentComplete = str(int(((self.steps/float(self.size))*1000))/10.0)
                res+=self.green+precentComplete+"%"+self.reset
                ex += len(self.green+self.reset)
            elif "fraction" in self.style:
                res+=":"
                res+=" "*max(1,30-len(res)+ex)
                res+=self.green+str(self.steps)+"/"+str(self.size)+self.reset
                ex += len(self.green+self.reset)
            elif "noProg" in self.style:
                res+=self.spinner[self.spinner_step]
            if self._estimate and not self._done:
                res+=" "*max(1,48-len(res)+ex)
                elapsed = time.time() - self.start_time
                estimate = elapsed/float(self.steps)*self.size - elapsed
                res+=self.yellow+_parseTime(estimate)+self.reset
                ex += len(self.yellow+self.reset)
            res = ShStr(res)
            res.ex = ex
            return res

    def subTracker(self,*args,**kwargs):
        with self.lock:
            kwargs["parent"] = self
            kwargs["lock"] = self.lock
            self.children.append(_Tracker(*args,**kwargs))
            return self.children[-1]

    def remove(self,subtracker):
        with self.lock:
            self.children.remove(subtracker)
            return self

    def base_tracker(self):
        base = self
        while base.parent!=None:
            base = base.parent
        return base

    def display(self):
        self.parent.display()
        return self

    def auto_display(self,rate=None):
        self.parent.auto_display(rate=rate)
        return self

    def out(self,content):
        self.parent.out(content)
        return self



class MultiTracker(_Tracker):
    def __init__(self, *args,**kwargs):
        super(MultiTracker, self).__init__(*args,**kwargs)
        self.printer = MultiPrinter()
        self.auto_rate = 0

    def auto_display(self,rate=None):
        next_time = time.time()+self.auto_rate
        with self.lock:
            if self._done:
                self.display()
            if rate != None: #Change rate if a rate has been set (0 means off)
                if rate==0:
                    self._auto_display = False
                else:
                    self.auto_rate = rate
                    next_time = time.time()+self.auto_rate
                    if self._auto_display:
                        return self
                    else:
                        self._auto_display = True
            if self._auto_display:
                self.display()
                durr = max(0,next_time-time.time())
                next_display = threading.Timer(durr,self.auto_display)
                next_display.daemon = True
                next_display.start()
            return self

    def freeze(self):
        with self.lock:
            print("\n"*(len(self.display(depth=1))-1))
            return self

    def done(self):
        with self.lock:
            super(MultiTracker, self).done()
            self.display()
            self.freeze()
            return self

    def display(self,to_disp=None,depth=0):
        lines = []
        if to_disp==None: to_disp = [self]
        for tracker in to_disp:
            tr = tracker.__str__()
            tString = ShStr("    "*depth+tr)
            tString.ex = tr.ex
            lines += [tString] + self.display(tracker.children,depth+1)
        if depth == 0:
            # adjusted_lines = []
            # offset = 0
            # for i,line in enumerate(lines):
            #     if i!=0:
            #         offset+=lines[i-1].ex
            #         adjusted_lines.append(" "*offset+line)
            #     else:
            #         adjusted_lines.append(line)
            self.printer.lines.change(lines)
        else:
            return lines
        return self

    def out(self,content):
        self.printer.out(content)
        return self

        

def _parseTime(seconds):
    if seconds < 60: return str(int(seconds))+"s"
    elif seconds == float('inf'): return "???"
    elif seconds < 3600: return str(int(seconds/60))+"m."+_parseTime(seconds%60)
    elif seconds < 86400: return str(int(seconds/3600))+"h."+_parseTime(seconds%3600)
    else: return str(int(seconds/86400))+"d."+_parseTime(seconds%86400)

class ShStr(str):
    """docstring for ClassName"""
    def __init__(self,*args,**kwargs):
        super(ShStr, self).__init__(*args,**kwargs)
        self.ex = 0
    def __len__(self):
        return len(str(self))-self.ex
    