% E-Flux 2 Script from publication was modified to include structure of iCC_CCBAU45436 and transcriptome standards 

%load a genome-scale metabolic model
model = readCbModel(strjoin({'..', 'Data', 'iCC541.mat'}, filesep));

% succinate production blocked
model = changeRxnBounds(model, {'EX_succ_e'}, 0, 'u');
% do loopless FVA to obtain tight bounds
if ~exist('fvaBounds.mat', 'file')
    modelFVA = model;
    % relax any required fluxes
    modelFVA.lb(modelFVA.lb > 0) = 0;
    modelFVA.ub(modelFVA.ub < 0) = 0;
    [minFlux, maxFlux] = fluxVariability(modelFVA, 0, 'max', [], 2, 0);
    assert(all(minFlux <= 0) & all(maxFlux >= 0))
    save('fvaBounds.mat', 'minFlux', 'maxFlux')
else
    d = load('fvaBounds.mat');
    [minFlux, maxFlux] = deal(d.minFlux, d.maxFlux);
end
%% Prepare for mapExpressionToReactions
fbamodel = model;
% I need to change LB and UB according to Eflux requirements: 
% LB = -10000 || 0 and  UB = 10000 || 0 will be used
for i = 1:length(fbamodel.rxns)
    rxnInd = findRxnIDs(fbamodel,fbamodel.rxns{i,1});
    
    if fbamodel.lb(rxnInd) < 0 
        fbamodel.lb(rxnInd) = -10000;
    end
end

fbamodel.ub(:) = 10000;    

if (~all(fbamodel.lb == 0 | fbamodel.lb == -10000) ...
    && ~all(fbamodel.ub == 0 | fbamodel.ub == 10000))
    disp('FBA model bounds must be 0 or +/- Infinity');
    return;
end

nrxn = length(fbamodel.rxns);
fbamodel.present = true(nrxn, 1); 

%% Load expression Data% W05: gene scores
load(strjoin({'..', 'Data', 'W05_geneValues.mat'}, filesep));

expressionData_W05.gene = fbamodel.genes;
expressionData_W05.value = W05_geneValues_v17;

% mapExpressionToReactions function determines the expression data associated to each reaction present in
% the model following standards

[expressionRxns_W05, parsedGPR_W05] = mapExpressionToReactions(fbamodel, expressionData_W05);

% C08: gene scores
load(strjoin({'..', 'Data', 'C08_geneValues.mat'}, filesep));

expressionData_C08.gene = fbamodel.genes;
expressionData_C08.value = C08_geneValues_v17;

[expressionRxns_C08, parsedGPR_C08] = mapExpressionToReactions(fbamodel, expressionData_C08);

%% E-flux using FVA bounds
% for each reaction, either both have expression values or none of them have expression values

[lbC, lbW] = deal(minFlux);
[ubC, ubW] = deal(maxFlux);

rxnWtExpression = expressionRxns_C08 >0 & expressionRxns_W05>0;

% for each reaction, get the condition with lower expression 
[~, ind] = min([expressionRxns_W05, expressionRxns_C08], [], 2);
% for each reaction, set the bounds for the condition with lower expression
% as bound * expression ratio
% reactions with lower expression in W05
rxnLowExprW05 = rxnWtExpression & ind == 1;
lbW(rxnLowExprW05) = minFlux(rxnLowExprW05) .* expressionRxns_W05(rxnLowExprW05) ./ expressionRxns_C08(rxnLowExprW05);
ubW(rxnLowExprW05) = maxFlux(rxnLowExprW05) .* expressionRxns_W05(rxnLowExprW05) ./ expressionRxns_C08(rxnLowExprW05);
% reactions with lower expression in C08
rxnLowExprC08 = rxnWtExpression & ind == 2;
lbC(rxnLowExprC08) = minFlux(rxnLowExprC08) .* expressionRxns_C08(rxnLowExprC08) ./ expressionRxns_W05(rxnLowExprC08);
ubC(rxnLowExprC08) = maxFlux(rxnLowExprC08) .* expressionRxns_C08(rxnLowExprC08) ./ expressionRxns_W05(rxnLowExprC08);

[modelW05, modelC08] = deal(model);
[modelW05.lb, modelW05.ub, modelC08.lb, modelC08.ub] = deal(lbW, ubW, lbC, ubC);
[solution1, solution2, totalFluxDiff] = optimizeTwoCbModels(modelW05, modelC08, 'max', false, true);
 
[f1, f2] = deal(struct());
for j = 1:numel(model.rxns)
    f1.(model.rxns{j}) = solution1.x(j);
    f2.(model.rxns{j}) = solution2.x(j);
end
f = fopen('fluxW05.json','w'); 
fprintf(f, jsonencode(f1));
fclose(f);
f = fopen('fluxC08.json','w'); 
fprintf(f, jsonencode(f2));
fclose(f);
