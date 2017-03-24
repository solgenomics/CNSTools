import time
import os, sys, threading

class MultiPrinter(object):
    """docstring for MultiPrinter"""
    def __init__(self, lines=None):
        self.lines = AccessAwareList(self,lines) if lines!=None else AccessAwareList(self,[])
        self.ms_size = 0 
        self.lock = threading.RLock()
        self.did_print = False

    def out(self,content):
        with self.lock:
            if content != "":
                orig = self.lines
                self.lines.change([" " for line in orig])
                self.print_members()

                ps = str(content)
                if ps[-2:]!="\n":
                    ps+="\n"
                sys.stdout.write(ps)

                self.lines = orig
                self.print_members()

    def print_members(self):
        with self.lock:
            if len(self.lines)==0:
                if self.did_print:
                    sys.stdout.write("\n\033[F\x1b[2K")
                    sys.stdout.flush()
                    self.did_print = False
                return
            width = get_terminal_size()[1]
            member_string = ""
            for line in self.lines:
                member_string += line
                padding = " "*width
                if len(line)>0:
                    overhang = len(line)%width
                    if overhang !=0:
                        padding = " "*(width-(len(line)%width))
                    else:
                        padding = ""
                member_string+=padding
            new_ms_size = len(member_string)//width-1
            sys.stdout.write("\r")
            if new_ms_size<self.ms_size:
                member_string+=" "*width*(self.ms_size-new_ms_size)
                sys.stdout.write(member_string+"\033[F"*self.ms_size)
            else:
                sys.stdout.write(member_string+"\033[F"*new_ms_size)
            self.ms_size = new_ms_size
            sys.stdout.flush()
            self.did_print = True

try:
    import fcntl, termios, struct
    def get_terminal_size(fd=1):
        """
        http://blog.taz.net.au/2012/04/09/getting-the-terminal-size-in-python/

        Returns height and width of current terminal. First tries to get
        size via termios.TIOCGWINSZ, then from environment. Defaults to 25
        lines x 80 columns if both methods fail.

        :param fd: file descriptor (default: 1=stdout)
        """
        try:
            return struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ, '1234'))
        except:  
            return (25, 80)

except:
    def get_terminal_size(fd=1):
        try:
            return (os.environ['LINES'], os.environ['COLUMNS'])
        except:  
            return (25, 80)


class AccessAwareList(list):
    def __init__(self, parent, obj):
        super(AccessAwareList, self).__init__(obj)
        self.parent = parent
        self.lock = threading.Lock()

    def __setitem__(self, key, value):
        with self.lock:
            super(AccessAwareList, self).__setitem__(key, value)
            self.parent.print_members()

    def change(self,lines):
        with self.lock:
            self.parent.lines = AccessAwareList(self.parent,lines)
            self.parent.print_members()
