To reproduce the figures of the article, proceed as follows.
The subfolder "reproducefigures" contains all the scripts.

The scripts "runthatfirst.ar1.generatedata.R" and "runthatfirst.ar5.generatedata.R"
generate data for the other scripts. That is, synthetic datasets from the hidden auto-regressive
models in dimensions 1 and 5 respectively.

The other scripts are named by figures.

# run the scripts 
figure1.likelihoodprofiles.run.d1.R
figure2-3.likelihoodprofiles.run.d5.R
# and then Figures 1-3 are produced by 
figure1-3.likelihoodprofiles.plots.R

# run the script
figure4.ar5.finitediff.plots.R
# and then Figure 4 is created by 
figure4.ar5.finitediff.run.R

# run the scripts
figure5.ar.kfmh.run.R
figure5.ar5T1000.pmmh.N64.R
figure5.ar5T1000.pmmh.N128.R
figure5.ar5T1000.pmmh.N256.R
# and then Figure 5 is produced by 
figure5.ar5T1000.acceptancerates.R
# and then Figure 6 is produced by
figure6.ar5T1000.posteriordensity.R


