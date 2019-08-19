The provided scripts require the following tools to be available:

Cobra Toolbox (https://github.com/opencobra/cobratoolbox), E-Flux algorithm ( from https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0157101 ),  Standardization algorithm ( from
https://academic.oup.com/bioinformatics/article/33/22/3603/3965323 )and Gurobi (Academic licenses are available free of charge via the academic program). The solver used in for this code was Gurobi v 8.0.1.


Function Descriptions:   
EfluxFVA.m - Generates lower and upper bounds using gene score data and run Eflux.
Yield.m - Calculate yields. Model needs to be normalized to calculate yields.
Validation.m - Model against published experimental data.
Sensitivity.m - Evolutionary simulation.



