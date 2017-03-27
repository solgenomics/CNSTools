import os, subprocess
from multitracker import MultiTracker

def call_command(command,shell=False,parent=None,tracker_name=None, env=os.environ,stdout=False,stderr=True):
    null_file = open(os.devnull, 'w')
    tracker = None  
    if tracker_name:
        if parent!=None:
            tracker = parent.subTracker(tracker_name,1,estimate=False,style="noProg")
        else:
            tracker = MultiTracker(tracker_name,1,estimate=False,style="noProg").auto_display(1)
        tracker.display()
    kwargs = {"env":env, "stdout":subprocess.PIPE if stdout else null_file, "stderr":subprocess.PIPE if stderr else null_file}
    if shell==True:
        proc = subprocess.Popen(" ".join(command),shell=True,**kwargs)
    else:
        proc = subprocess.Popen(command,**kwargs)
    proc.wait()
    output = print_cmd_err(proc,stdout,stderr)
    if tracker:
        tracker.base_tracker().printer.out(output)
        tracker.done().display()
    else:
        print(output)        
    null_file.close()
    return proc

def call_commands_async(command_list,num,shell=False,parent=None,tracker_name=None, env=os.environ,stdout=False,stderr=True,estimate=False):
    null_file = open(os.devnull, 'w')
    process_list = []
    finished = []
    if len(command_list)<num: num = len(command_list)
    tracker = None  
    if tracker_name:
        if parent!=None:
            tracker = parent.subTracker(tracker_name,len(command_list),estimate=estimate,style="fraction")
        else:
            tracker = MultiTracker(tracker_name,len(command_list),estimate=estimate,style="fraction").auto_display(1)
        tracker.display()   
        tracker.toggle()
    for command in command_list:
        #print " ".join(command)
        kwargs = {"env":env, "stdout":subprocess.PIPE if stdout else null_file, "stderr":subprocess.PIPE if stderr else null_file}
        if shell==True:
            process_list.append(subprocess.Popen(" ".join(command),shell=True,**kwargs))
        else:
            process_list.append(subprocess.Popen(command,**kwargs))
        while len(process_list) >= num:
            pid,exitstat = os.waitpid(-1,0)
            for i in range(len(process_list)-1,-1,-1):
                if process_list[i].pid == pid or process_list[i].poll() != None:
                    finished.append(process_list.pop(i))
                    output = print_cmd_err(finished[-1],stdout,stderr)
                    if tracker: 
                        tracker.base_tracker().printer.out(output)
                        tracker.step()
                    else:
                        print(output)
    while len(process_list) > 0:
        pid,exitstat = os.waitpid(-1,0)
        for i in range(len(process_list)-1,-1,-1):
            if process_list[i].pid == pid or process_list[i].poll() != None:
                finished.append(process_list.pop(i))
                output = print_cmd_err(finished[-1],stdout,stderr)
                if tracker: 
                    tracker.base_tracker().printer.out(output)
                    tracker.step()
                else:
                    print(output)
    if tracker: 
        tracker.done().display()
    null_file.close()
    return finished

def print_cmd_err(p,stdout,stderr):
    stdout_text,stderr_text = p.communicate()
    ps = ""
    if stdout and stdout_text:
        ps += "\x1b[34m%s\x1b[0m"%stdout_text.strip() #prints stdout in blue
    if stderr and stderr_text:
        if ps!="":
            ps+="\n"
        ps += "\x1b[31m%s\x1b[0m"%stderr_text.strip() #prints stderr in red
    return ps