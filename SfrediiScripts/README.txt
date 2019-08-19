# Scripts provided for the Paper "Metabolic analyses of nitrogen fixation in the soybean micro-symbiont Sinorhizobium fredii using constraint-based modeling"
By Carolina A. Contador, Siu-Kit Lo, Siu H. J. Chan, and Hon-Ming Lam 

## The provided scripts require the following tools to be available:

Cobra Toolbox ([available here](https://github.com/opencobra/cobratoolbox)) and Gurobi (Academic licenses are available free of charge via the ([Gurobi academic program](https://www.gurobi.com/academia/academic-program-and-licenses/)). The solver used in the paper was Gurobi v8.1.0.  


## Matlab Function Descriptions:
All Matlab Scripts are contained in the *Matlab* folder.
EfluxFVA.m - Generates lower and upper bounds using gene score data and run Eflux.
Yield.m - Calculate yields. Model needs to be normalized to calculate yields.
Validation.m - Model against published experimental data.
Sensitivity.m - Evolutionary simulation.


## Memote Descriptions:
All Files required to run memote are contained in the *Memote* folder.
Memote ([available here](https://github.com/opencobra/memote)) 
## Escher Descriptions:
All Files required to create the flux distributions using Escher are contained in the *Escher* folder.
Escher ([available here](https://escher.github.io/)) 
## S. Meliloti models:
xml files were downloaded from supplementary material of ([iHZ565](https://doi.org/10.1371/journal.pone.0031287.s002)) and ([iGD1575] (https://www.nature.com/articles/ncomms12219)).

## Please cite:
Contador, C.A., Lo, S-K., Chan, S.H.J., Lam, H-M. (2019). Metabolic analyses of nitrogen fixation in the soybean micro-symbiont Sinorhizobium fredii using constraint-based modeling. (submitted).
