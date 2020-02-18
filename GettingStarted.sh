# This shows how to structure inputs for GenoML CLI 

## SETTING UP A VIRTUAL ENVIRONMENT 
# Making a virtual environment
conda create -n GenoML python=3.7

# Activating and changing directories to environment
conda activate GenoML
    # Deactivating a conda environment 
        # conda deactivate ENV_NAME 
    # Removing a conda environment 
        # conda env remove -n ENV_NAME

# Installing from a requirements file using pip
pip install -r requirements.txt

# Install the package at this path
pip install .

# Saving out environment requirements to a .txt file
#pip freeze > requirements.txt

# Removing a conda virtualenv
#conda remove --name GenoML --all

# Run GenoML 
GenoML
#usage: GenoML [-h] [--prefix PREFIX] [--geno GENO] [--addit ADDIT]
            #  [--pheno PHENO] [--gwas GWAS] [--p P] [--vif VIF] [--iter ITER]
            #  [--impute IMPUTE] [--rank_features {skip,run}]
            #  [--max_tune MAX_TUNE] [--n_cv N_CV]
            #  {discrete,continuous} {supervised,unsupervised} {train,tune}
#GenoML: error: the following arguments are required: data, method, mode

# Running the munging script
GenoMLMunging --prefix outputs/test_discrete_geno \
--geno examples/training \
--pheno examples/training_pheno.csv 

# Running the munging script with VIF
GenoMLMunging --prefix outputs/test_discrete_geno \
--geno examples/training \
--pheno examples/training_pheno.csv \
--vif 5 \
--iter 1

# Running the discrete supervised training script
GenoML discrete supervised train \
--prefix outputs/test_discrete_geno \
--rank_features run

# Running the discrete supervised tuning script


# Running the continuous supervised training script
GenoML continuous supervised train \
--prefix outputs/test_discrete_geno \
--rank_features run

# Running the continuous supervised tuning script 

