#David Lyon    
                                                                                                                                                                                 
import time

class Progress_tracker():
    """docstring for Progress"""
    def __init__(self,name,size,startIt=False):
        self.name = name
        self.size = size
        self.prevStep = None
        self.timeElapsed = 0
        self.steps = 0
        self.prog = 0.000001
        self.prevDisplay = 0
        self.autoDisplay = False
        self.displayRate = 0
        if startIt: self.start()
    def start(self):
        self.prevStep = time.time()
    def step(self,size=1):
        if not self.prevStep: self.start()
        elif self.steps<self.size:
            thisStep = time.time()
            self.timeElapsed+= thisStep-self.prevStep
            self.prevStep = thisStep
            self.steps += size
            if self.steps>self.size: self.steps=self.size
            self.prog = (self.steps/float(self.size))
            if thisStep-self.prevDisplay>self.displayRate:
                self.display()
    def display(self, estimate=True, rate=None):
        if rate>0:
            self.autoDisplay = True
            self.displayRate = rate
        elif rate:
            self.autoDisplay = False
            self.displayRate = rate
        out = self.name+': '
        if self.steps<self.size:
            out+=str(int((self.prog*100)))+"%"
            if estimate:
                out+=" ("+self.parseTime(int(self.timeElapsed/self.prog*(1-self.prog)))+" remaining)"
            else:
                out+="..."
        else:
            out+='Finished'
        self.prevDisplay = time.time()
        print out
    def parseTime(self,seconds):
        if seconds < 60: return str(seconds//1)+"s"
        elif seconds < 3600: return str(seconds//60)+"m."+self.parseTime(seconds%60)
        elif seconds < 86400: return str(seconds//3600)+"h."+self.parseTime(seconds%3600)
        else: return str(seconds//86400)+"d."+self.parseTime(seconds%86400)




