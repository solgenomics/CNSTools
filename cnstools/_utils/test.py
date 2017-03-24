import curses
import time
import multitracker
import random

t = multitracker.MultiTracker("Testing whooo",10000,style="percent").auto_display(1);
time.sleep(0.5)
# t.subTracker("Testing 1",1,style="fraction",estimate=False)
# k = t.subTracker("Testing 2",1,style="silent",estimate=False)
# t.subTracker("Testing 3",1,style="silent",estimate=False)
# t.subTracker("Testing 4",1,style="silent",estimate=False)
i = 0
while not t._done:
    i+=1
    time.sleep(0.5)
    t.step(random.randint(1,1000))
t.remove(k)
time.sleep(10)
t.freeze()

# def report_progress(filename, progress):
#     """progress: 0-10"""
#     stdscr.addstr(0, 0, "Moving file: {0}".format(filename))
#     stdscr.addstr(1, 0, "Total progress: [{1:10}] {0}%".format(progress * 10, "#" * progress))
#     print("Hi")
#     stdscr.refresh()

# if __name__ == "__main__":
#     stdscr = curses.initscr()
#     curses.noecho()
#     curses.cbreak()

#     try:
#         for i in range(10):
#             report_progress("file_{0}.txt".format(i), i+1)
#             time.sleep(0.5)
#     finally:
#         curses.echo()
#         curses.nocbreak()
#         curses.endwin()