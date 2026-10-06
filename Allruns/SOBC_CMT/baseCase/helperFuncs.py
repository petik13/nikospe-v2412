import numpy as np
import os
import subprocess
import logging
from typing import List
import shutil

def format_time(t: float) -> str:
    return f"{t:.6f}"

class SimpleLogger:
    """A very simple logger that definitely works"""
    def __init__(self, log_file: str = 'log.txt'):
        self.log_file = os.path.abspath(log_file)
        print(f"SimpleLogger will write to: {self.log_file}")
        # Create/test the file immediately
        with open(self.log_file, 'a') as f:
            f.write("=== Logger started ===\n")
            f.flush()
        print(f"Log file created successfully: {self.log_file}")
    
    def _write(self, level: str, message: str):
        line = f"[{level}] {message}\n"
        print(line.strip())  # Also print to console
        with open(self.log_file, 'a') as f:
            f.write(line)
            f.flush()
    
    def info(self, message: str):
        self._write('INFO', message)
    
    def warning(self, message: str):
        self._write('WARNING', message)
    
    def error(self, message: str):
        self._write('ERROR', message)

def setup_logger(log_file: str = 'log.txt'):
    """Setup and return a simple, reliable logger"""
    return SimpleLogger(log_file)

def update_control_dict(start: float, end: float, *, logger):
    path = os.path.join('system', 'controlDict')
    if not os.path.isfile(path):
        logger.warning(f"controlDict not found at {path}, skipping update.")
        return
    # in-place edit minimal (sed-like) using Python
    with open(path, 'r') as f:
        lines = f.readlines()
    def repl(line: str, key: str, value: float) -> str:
        if line.strip().startswith(key):
            return f"{key}\t{value};\n"
        return line
    new_lines = [repl(repl(ln, 'startTime', start), 'endTime', end) for ln in lines]
    with open(path, 'w') as f:
        f.writelines(new_lines)

def console(msg: str):
    print(msg)

def run_cmd(cmd: List[str], *, logger, check=True, env=None) -> int:
    """Run a command streaming stdout/stderr to the logfile in real-time."""
    shown = ' '.join(cmd)
    logger.info(f"$ {shown}")
    
    try:
        # Start the process with pipes for real-time output capture
        process = subprocess.Popen(
            cmd, 
            stdout=subprocess.PIPE, 
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1,
            env=env
        )
        
        # Stream output line by line to logger
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                # Remove trailing newline since logger adds its own
                line = output.rstrip('\n\r')
                if line.strip():  # Only log non-empty lines
                    # Write directly to log file without [INFO] prefix for command output
                    with open(logger.log_file, 'a') as f:
                        f.write(f"{line}\n")
                        f.flush()
        
        # Get the return code
        return_code = process.poll()
        
    except KeyboardInterrupt:
        logger.warning("Command interrupted by user (KeyboardInterrupt)")
        process.terminate()
        raise
    except Exception as e:
        logger.error(f"Error running command: {e}")
        return 1
        
    if check and return_code != 0:
        logger.error(f"Command failed (code {return_code}): {shown}")
        raise SystemExit(return_code)
    
    return return_code

def list_time_dirs(case_root: str='.') -> List[str]:
    dirs = []
    for name in os.listdir(case_root):
        if name in ('0', '0.org', 'constant', 'system', 'postResults'): continue
        p = os.path.join(case_root, name)
        if os.path.isdir(p):
            try:
                float(name)
            except ValueError:
                continue
            dirs.append(name)
    return sorted(dirs, key=lambda s: float(s))

def stash_results(at_time: float):
    label = format_time(at_time)
    dest = os.path.join('postResults', label)
    os.makedirs(dest, exist_ok=True)
    # Move ALL numeric time directories (excluding 0, 0.org, constant, system, postResults)
    moved = 0
    for tdir in list_time_dirs('.'):
        shutil.move(tdir, dest)
        moved += 1
    # Copy 0, constant, system snapshots
    for d in ('0', 'constant', 'system'):
        if os.path.isdir(d):
            target = os.path.join(dest, d)
            if os.path.exists(target):
                shutil.rmtree(target)
            if d == '0':
                shutil.copytree('0', target, symlinks=False)
            else:
                shutil.copytree(d, target, symlinks=False)
    open(os.path.join(dest, 'foam.foam'), 'a').close()