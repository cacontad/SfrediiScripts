load('nodules_models.mat')
model = modelW05;

model_goal = modelC08;

objFlux = [];
target = [];
model2 = model;


for i = 1:10000
    % pick random reaction
    metrxn = model.rxns{randi(numel(model.rxns), 1)};
    rxnInd = findRxnIDs(model,metrxn);
    % bounds
    W05_lb = model2.lb(rxnInd);
    W05_ub = model2.ub(rxnInd);
    C08_lb = model_goal.lb(rxnInd);
    C08_ub = model_goal.ub(rxnInd);
    lb1 = min(W05_lb,C08_lb);
    lb2 = max(W05_lb,C08_lb);
    ub1 = min(W05_ub,C08_ub);
    ub2 = max(W05_ub,C08_ub);
    lb = lb1 + (lb2-lb1).*rand(1,1);
    ub = ub1 + (ub2-ub1).*rand(1,1);
    model2.lb(rxnInd) = lb;
    model2.ub(rxnInd) = ub;
    s = optimizeCbModel(model2);
    objFlux(i) = s.f;
    target(i) = rxnInd;
end


objFlux = objFlux';

plot(objFlux,'*')
xlabel('Iteration');
ylabel('Symbiosis flux (mmol gDW^{-1} hr^{-1})');
