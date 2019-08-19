load('iCC541.mat')
model = iCC541;


FBAsolution_wildtype = optimizeCbModel(model,[],'one');


% Validation script for iCC541 reconstruction
% Scenarios tested were extracted from published literature


%% check phb production
[modelphb,rxnNames] = addDemandReaction(model,'phb[c]');
modelphb = changeObjective(modelphb,'DM_phb[c]',1);
FBAsolution_phb = optimizeCbModel(model,[],'one');


%% check pyrC deletion: Fix- phenotype
% AB395_0000120 = AB395_RS00565
geneID = findGeneIDs(model,'AB395_0000120');
geneList = model.genes(geneID,1);
[model_mutPyrC] = deleteModelGenes(model,geneList);
FBAsolution_pyrC = optimizeCbModel(model_mutPyrC,[],'one');

%% check pyrF deletion: Fix- phenotype
% AB395_00004141 = AB395_RS20020
geneID = findGeneIDs(model,'AB395_00004141');
geneList = model.genes(geneID,1);
[model_mutPyrF] = deleteModelGenes(model,geneList);
FBAsolution_pyrF = optimizeCbModel(model_mutPyrF,[],'one');

%% C4 requirement
model_C4 = model;
% Succinate_uptake = -1.38
model_C4 =changeRxnBounds(model_C4,'EX_succinate',0,'l');
% Malate uptake = -1.44
model_C4 =changeRxnBounds(model_C4,'EX_malL',0,'l');
model_C4 =changeRxnBounds(model_C4,'EX_gluL',-0.6,'l');
FBAsolution_C4 = optimizeCbModel(model_C4,[],'one');

%% Check symbiosis sensitivity on malate and succinate
[sucUt, malUt] = deal(0:0.05:1, 0:0.04:1); 
objSucMal = zeros(numel(sucUt), numel(malUt));
for i = 1:numel(sucUt)
    for j = 1:numel(malUt)
        model2 = changeRxnBounds(model, {'EX_malL', 'EX_succinate'}, ...
            -[malUt(j); sucUt(i)], 'l');
        sIJ = optimizeCbModel(model2);
        objSucMal(i, j) = sIJ.obj;
    end
end
subplot(2,2,3);
surf(sucUt, malUt, objSucMal');
xlabel('succinate uptake (mmol gDW^{-1} hr^{-1})');
ylabel('malate uptake (mmol gDW^{-1} hr^{-1})');
zlabel('Symbiosis flux (mmol gDW^{-1} hr^{-1})');

%% check purL deletion: Fix- phenotype
% AB395_00001672 = AB395_RS08020
geneID = findGeneIDs(model,'AB395_00001672');
geneList = model.genes(geneID,1);
[model_mutPurL] = deleteModelGenes(model,geneList);
FBAsolution_purL = optimizeCbModel(model_mutPurL,[],'one');

%% check purQ deletion: Fix- phenotype
% AB395_00001668 = AB395_RS08000
geneID = findGeneIDs(model,'AB395_00001668');
geneList = model.genes(geneID,1);
[model_mutPurQ] = deleteModelGenes(model,geneList);
FBAsolution_purQ = optimizeCbModel(model_mutPurQ,[],'one');

%% check idhA: inositol dehydrogenase
% AB395_00002480 = AB395_RS11960
geneID = findGeneIDs(model,'AB395_00002480');
geneList = model.genes(geneID,1);
[model_mutidh] = deleteModelGenes(model,geneList);
model_mutidh = changeRxnBounds(model_mutidh,'EX_inost',0,'l');
FBAsolution_idh = optimizeCbModel(model_mutidh,'max','one');

%% check znuA deletion: Zinc ABC transport (high-affinity transport system)

%r_0522: zinc diffusion transport rxn
modelZn = changeRxnBounds(model,'r_0522',-0.0001,'b');

%ZnuA: AB395_RS08730
% AB395_00001824 = AB395_RS08730
geneID = findGeneIDs(model,'AB395_00001824');
geneList = model.genes(geneID,1);
[model_mutZn] = deleteModelGenes(modelZn,geneList);
FBAsolution_Zn = optimizeCbModel(model_mutZn,[],'one');

%% check ptsSCAB mutant: high-affinity transporter
% in-frame deletion mutants of ptsSCAB
%r_0546: phosphate transport rxn (low-affinity)
modelP = model;
modelP = changeRxnBounds(modelP,'r_0546',-0.001,'l');
modelP = changeRxnBounds(modelP,'r_0546',0.001,'u');
modelP = changeRxnBounds(modelP,'r_0313',0,'b');
FBAsolution_P = optimizeCbModel(modelP,[],'one');

%% check cobO mutant: cobalamin

%modelCob = model;
%modelCob = changeRxnBounds(modelCob,'EX_cbl1',-1000,'l');
% modelCob = changeRxnBounds(modelCob,'r_0254',0);
% modelCob = changeRxnBounds(modelCob,'r_0035',0);
% cobO: AB395_RS08990
% AB395_00001876 = AB395_RS08990
geneID = findGeneIDs(model,'AB395_00001876');
geneList = model.genes(geneID,1);
[model_mutCob] = deleteModelGenes(model,geneList);
FBAsolution_Cob = optimizeCbModel(model_mutCob,'max','one');

