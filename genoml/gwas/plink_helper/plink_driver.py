import os
import glob
import subprocess

# local import
from plink_helper.plink_logger import Logger

class Driver(Logger):
    def __init__(self, geno_path):
        super().__init__(geno_path)
        
    def run_cmds(self, cmds_list, step):
        log = self.logging()
        log.write(step + " WITH THE FOLLOWING COMMANDS:")
        log.write("\n")
        log.write("\n")
        
        
        for cmd in cmds_list:
            log.write(cmd)
            log.write("\n")
            log.write("\n")
            
            # check length of list of files ending in ".log"
            logs_count = len(sorted(glob.glob(self.out_path + '*.log'), key=os.path.getmtime))
            subprocess.run(cmd, shell=True)
            
            # check length of new list of files ending in ".log". if longer than original, append newest .log file to running logfile
            # this indicates new log created
            new_logs_count = len(sorted(glob.glob(self.out_path + '*.log'), key=os.path.getmtime))
            
            if new_logs_count > logs_count:
                cmd_log = sorted(glob.glob(self.out_path + '*.log'), key=os.path.getmtime)[-1]
                
                new_log = open(cmd_log, "r")
                new_log_read = new_log.read()
                new_log.close()
                
                log.write(new_log_read)
                log.write("\n")
                log.write("***********************************************")
                log.write("\n")
                log.write("\n")

        log.close()
        