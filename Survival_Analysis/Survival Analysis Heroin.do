*```````````````````````````````````````````````````````````````````````````````
*                           SURVIVAL ANALYSIS
*```````````````````````````````````````````````````````````````````````````````
cd "D:\Field Epidemiology 2026\Introduction to Biostatistics II"
clear


*             Description
desc 
codebook

*             Prepare the data

label define binary 0 "No" 1 "Yes"

*            Survival Rendering
stset time, failure(status)
