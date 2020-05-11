import subprocess
import argparse
import os
import requests
import json
import time
import glob

#local imports
from plink_helper.plink_driver import Driver
import genoml.dependencies

# parser = argparse.ArgumentParser(description='Arguments for Genotyping Imputation (data in Plink .bim/.bam/.fam format)')
# parser.add_argument('--geno', type=str, default='nope', help='Genotype: (string file path). Path to PLINK format genotype file, everything before the *.bed/bim/fam [default: nope].')
# parser.add_argument('--out', type=str, default='out', help='Prefix for output (including path)')

# args = parser.parse_args()

# geno = args.geno
# out = args.out
geno = '/data/vitaled2/test_data/PDBP/PDBP_het_call_rate_sex_relatedness_variant_final'
out = '/data/vitaled2/test_data/PDBP/'



################################################

class Impute(Driver):
    def __init__(self, geno_path):
        super().__init__(geno_path)
        self.vcf_list_for_impute = [geno_path + "_pre_impute" + "_chr" + str(i) + ".vcf.gz" for i in range(1,24)]
        self.plink_exec = genoml.dependencies.check_plink()
        
    def impute_prep_data(self):
        geno_path = self.geno_path
        out_path = self.out_path
        plink_exec = self.plink_exec
        step = "PREP PLINK FILES FOR IMPUTATION"
        os.chdir(out_path)
        # download file to check
        bash1 = "wget http://www.well.ox.ac.uk/~wrayner/tools/HRC-1000G-check-bim.v4.2.5.zip -P " + out_path

        bash2 = "unzip " + out_path + "HRC-1000G-check-bim.v4.2.5.zip -d " + out_path

        bash3 = "wget ftp://ngs.sanger.ac.uk/production/hrc/HRC.r1-1/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz -P " + out_path
        bash4 = "gunzip " + out_path + "HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz"
        # make your .frq file
        bash5 = f"{plink_exec} --bfile " + geno_path + " --freq --out " + geno_path

        bash6 = "perl " + out_path + "HRC-1000G-check-bim.pl -b " + geno_path + ".bim -f " + geno_path + ".frq -r " + out_path + "HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h"

        # then run to fix your data
        bash7 = "sh " + out_path + "Run-plink.sh"

        cmds = [bash1, bash2, bash3, bash4, bash5, bash6, bash7]

        self.run_cmds(cmds, step)
            

    def impute_make_vcf(self):
        geno_path = self.geno_path
        out_path = self.out_path
        plink_exec = self.plink_exec
        # then make vcf files
        step1 = "RECODE PLINK FILES TO VCF"
        
        cmds1 = [f"{plink_exec} --bfile " + geno_path + "-updated-chr" + str(i) + " --recode vcf --chr " + str(i) + " --out " + geno_path + "_chr" + str(i) for i in range(1,24)]   

        self.run_cmds(cmds1, step1)
        
        ## then sort and zip
        step2 =  "vcf-sort AND bgzip VCFS"
        cmds2 = ["vcf-sort " + geno_path + "_chr" + str(i) + ".vcf | bgzip -c > " + geno_path + "_pre_impute" + "_chr" + str(i) + ".vcf.gz" for i in range(1,24)]

        self.run_cmds(cmds2, step2)

        
    def check_impute_status(self, _key, _id):
        
        # imputation server url
        url = 'https://imputationserver.sph.umich.edu/api/v2'

        # add token to header (see authentication)
        headers = {'X-Auth-Token' : _key }

        # get all jobs
        r = requests.get(url + "/jobs", headers=headers)
        if r.status_code != 200:
            raise Exception('GET /jobs/ {}'.format(r.status_code))
        
        status = r.json()
        for stat in status['data']:
            if stat['id'] == _id:
                if stat['state'] == 1:
                    print("Launching Job:", stat['id'])
                elif stat['state'] == 2:
                    print("Running Job:", stat['id'])
                elif stat['state'] == 3:
                    print(stat['id'], "returned state '3', have a look at jobs on the web front for more information")
                elif stat['state'] == 5:
                    print(stat['id'], "has failed. consult docs on data input to ensure your vcfs are correct")
                elif stat['state'] == 4:
                    print(stat['id'], "COMPLETED!")
                
                return stat['state']
            
            else:
                pass
        

    def pull_imputed_data(self, _key, _id, pw):
        out_path = self.out_path
        geno_path = self.geno_path
        
        #create path for impute ouput files
        imputed_pathname = out_path + 'imputed'
        os.mkdir(imputed_pathname)
        os.chdir(imputed_pathname)
        
        # imputation server url
        url = 'https://imputationserver.sph.umich.edu/api/v2'

        # add token to header (see authentication)
        headers = {'X-Auth-Token' : _key }
        
        r = requests.get(url + "/jobs/" + _id, headers=headers)
        if r.status_code != 200:
            raise Exception('GET /jobs/ {}'.format(r.status_code))

        output_json = r.json()

        hashes_dict = {output_json['outputParams'][i]['id'] : output_json['outputParams'][i]['hash'] for i in range(len(output_json['outputParams']))}
        
        # run a curl for each
        curls = ['curl -sL https://imputationserver.sph.umich.edu/get/' + str(key) + '/' + str(hashes_dict[key]) + ' | bash' for key in hashes_dict]
        
        for curl in curls:
            print("Curling output data with the following command: " + curl)
            subprocess.run(curl, shell=True)
        print() 
        print("Finished Pulling Imputed Data!")
        print()
        
        #now unzip all ".zip" files (one for each chromosome)
        zip_list = glob.glob(out_path + 'imputed/*.zip')
        unzip_cmds = ['unzip -P ' + pw + ' ' + file for file in zip_list]
        
        for cmd in unzip_cmds:
            print("Unzipping: " + cmd)
            subprocess.run(cmd, shell=True)
        print("Finished Unzipping")
   

    def impute(self, key, input_population='eur', pw='imputer', vcf_list=None):
        geno_path = self.geno_path
        vcf_list = self.vcf_list_for_impute
        
        # imputation server url
        url = 'https://imputationserver.sph.umich.edu/api/v2'

        # add token to header (see Authentication)
        headers = {'X-Auth-Token' : key}

        open_vcfs = [open(vcf, 'rb') for vcf in vcf_list]
        
        files = set([('input-files-upload', vcf) for vcf in open_vcfs])

        data = {'input-mode' : 'imputation',
                'input-files-source': 'file-upload',
                'input-password': pw,
                'input-refpanel': 'apps@hrc-r1.1',
                'input-phasing': 'eagle',
                'input-population': input_population}

        r = requests.post(url + "/jobs/submit/minimac4", files=files, headers=headers, data=data)
        if r.status_code != 200:
            raise Exception('POST /jobs/submit/minimac4 {}'.format(r.status_code))
        
        impute_id = r.json()['id']
        message = r.json()['message']

        print(message)
        print(impute_id)
        print('***************************')
        print('* * * * * * * * * * * * * *') 
        
        imp_state = 0
        while imp_state < 3:
            time.sleep(60)
            os.system('clear')
            imp_state = self.check_impute_status(key, impute_id)
            
            if imp_state == 4:
                print("Pulling Completed Data from Imputation Server!")
                self.pull_imputed_data(key, impute_id, pw)
                
                
    def cleanup(self):
        geno_path = self.geno_path
        out_path = self.out_path
        print("CLEANING UP THE DIRECTORY OF INTERMEDIATE FILES")
        print("***********************************************")
        print()
        
        # make list of hrc files but keep LOG
        split_geno_files = glob.glob(geno_path + '*chr*')
        hrc_files = list(set(glob.glob(out_path + '*HRC*')) - set(glob.glob(out_path + 'LOG*HRC*')))
        updated_files = glob.glob(out_path + '*-updated.*')

        rm_files = split_geno_files + hrc_files + updated_files + [out_path + 'LICENSE.txt', out_path + 'Run-plink.sh']
        
        for file in rm_files:
            os.remove(file)
            

imputer = Impute(geno)
# imputer.impute_prep_data()
# imputer.impute_make_vcf()
imputer.impute(key=<put key here!!!!>)
