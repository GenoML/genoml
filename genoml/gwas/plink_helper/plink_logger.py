import os

class Logger:
    def __init__(self, geno_path):
        self.geno_path = geno_path
        self.out_path = self.geno_path.rsplit('/', 1)[0] + "/"
        print("PROCESSING THE FOLLOWING GENOTYPES:", geno_path)
        print()
        
        
    def logging(self):
        log_out_name = self.geno_path + ".PLINK_STEPS.log"
        print("logging to " + log_out_name)
        print("***********************************************")
        print()
        log = open(log_out_name, "a", newline='\n')
        log.write("***********************************************")
        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
        log.write("RUNNING PLINK COMMANDS FOR " + self.geno_path)
        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
        log.write("***********************************************")
        log.write("\n")
        log.write("\n")        
            
        return log
    
    def rm_log(self):
        log_out_name = self.geno_path + ".PLINK_STEPS.log"
        if os.path.exists(log_out_name):
            os.remove(log_out_name)
        