import subprocess
import multiprocessing
import time
import os

def process_monotor(delay):
    process = subprocess.Popen('ps -ef| grep run_regressions.R | grep -v ? | grep -v "grep" | wc -l', stdout=subprocess.PIPE,shell=True)
    p_count = int(list(process.communicate())[0].decode().rstrip())
    time.sleep(delay)
    return p_count

def load_process(cmd):
    os.system(cmd + '&')
    return None

def driver(commands,scaling_factor=1, delay=0):
    running_processes = 0
    for c in commands: 
        if running_processes <= 63: 
            running_processes = process_monotor(delay=delay)
            load_process(c)
            print(running_processes)
        while running_processes > 63:
            running_processes = process_monotor(delay=5)
    return None

if __name__ == '__main__':
    commands = []
    with open('commands') as f:
        for line in f:
            commands.append(line.rstrip())
    driver(commands)

