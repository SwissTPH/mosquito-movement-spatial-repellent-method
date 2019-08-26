# mosquito-movement-spatial-repellent-method
method for secondary analysis of trial data

Contains code from the paper <br>
<i>Can trials of spatial repellents be used to estimate mosquito movement?</i> Malinga J, Maia M, Moore S & Ross A. <br>
<i>Parasites & Vectors</i> 2019, in press

# cpp codes _ simulation

Contains the C++ code to perform the analysis on the trial data (which can be real or simulated).
There are separate files for the trials with/out additional seasonality data. 

Each input dataset has to be renamed to "aggTrapSimul.txt" for the code to run (aggTrapSimul.txt contains an example simulated dataset, which can be made using the R scripts in /r_codes). <br>

The output file results.txt gives the estimated baseline distribution of mosquitoes in the village (in the case of zero coverage of spatial repeelents) as a proportion for each house in the village. The parameter values for repFactor (the proportion of mosquitos that are diverted by spatial repellents of those in houses that use repellents), hseFactor (the proportion of mosquitoes diverted elsewhere as opposed to another house), lambda (the mean distance between households moved by the mosquitoes) and the log likelihood are given underneath.  

The fitting method uses the Nelder-Mead algorithm from the AS047 library by John Burkhardt. This came from:<br>
https://people.sc.fsu.edu/~jburkardt/f77_src/asa047/asa047.html


# r codes _ simulation

Contains the R scripts for simulating the trial scenarios (with/out additional seasonality data). The simulated datasets are used to test how well the analysis method would work for a trial with specific parameter values.

The parameter values and file names for the simulations can be edited at the end of each R script.
These R scripts can simulate datasets with different parameter combinations.
