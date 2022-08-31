# heteroplasmySEVar-and-Kimura-fits

Scripts associated with the manuscript XXX (2022). 

Script to plots the standard errors of the variance.

Some examples of how to call the script with arguments are those:
Rscript SEVarPlots.R 'real' 'Freyer' 'NA' -1 -1 1
Rscript SEVarPlots.R 'real' 'orgDatBroz.txt' 'NA' -1 -1 1
Rscript SEVarPlots.R 'synthetic' 'NA' 'rnorm' 0.5 0.1 3
Rscript SEVarPlots.R 'synthetic' 'NA' 'kimura' 0.6 0.9 3

in general, the input arguments have the below form
Rscript SEVarPlots.R    ['synthetic'|'real']       ['NA'|'HB'| 'LE'| 'Freyer'| 'Broz'| a custom path to a txt]        ['rbeta' | 'rlnorm' | 'kimura' | 'rnorm' | 'NA']        [distr parameter 1]        [distr parameter 2]   [# rows in the plot grid]

And scripts to calculate the p-vals.

Script pvalsCode.R is used to produce the p vals for the Freyer data. 
Rscript pvalsCode.R  'Offsprings' 100 10
where the first input argument can be PGC, Oocytes, or Offsprings (read from the attached txt files)the second is the number of MCMC samples and the third one is the number of repetitions per sample

Script pvalsCodeBinned.R is used for the binned data from Wonnapinij.
E.g., Rscript  pvalsCodeBinned.R  100 10 10
(first two arguments the same as above, the third argument is the number of resamples from the binned data)

Both of the p-vals scripts can also be executed in parallel with slight modifications (follow comments mase in the scripts).
