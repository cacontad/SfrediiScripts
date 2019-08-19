% load your model version
load('iCC541.mat');
model = iCC541;


% Symbiosis product molecular weight was set to 1 g/mmol using 
% Standardization algorithm ( from https://academic.oup.com/bioinformatics/article/33/22/3603/3965323 )

s = optimizeCbModel(model, 'max', 'one');
fprintf('Normalized symbiosis flux: %.4f\n', s.f);

rxnEx = find(sum(model.S ~= 0, 1) == 1);
[metEx, ~] = find(model.S(:, rxnEx));
[Masses, knownMasses, unknownElements, Ematrix, elements] = getMolecularMass(model.metFormulas(metEx), 0, true);
eleC = strcmp(elements, 'C');
utC = s.x(rxnEx) < 0 & Ematrix(:, eleC) > 0;
yieldC = s.f / abs(s.x(rxnEx(utC))' * Ematrix(utC, eleC)); 
fprintf('Symbiosis flux: %.4f g/C-mol\n', yieldC)


eleC = strcmp(elements, 'N');
utC = s.x(rxnEx) < 0 & Ematrix(:, eleC) > 0;
yieldN = s.f / abs(s.x(rxnEx(utC))' * Ematrix(utC, eleC)); 
fprintf('Symbiosis flux: %.4f g/N-mol\n', yieldN)
ratio = yieldC/yieldN;
fprintf('C/N: %.4f g/C-mol\n', ratio)


% W05
fprintf('Normalized symbiosis flux W05: %.4f\n', solution1.f);
rxnEx = find(sum(model.S ~= 0, 1) == 1);
[metEx, ~] = find(model.S(:, rxnEx));
[Masses, knownMasses, unknownElements, Ematrix, elements] = getMolecularMass(model.metFormulas(metEx), 0, true);
eleC = strcmp(elements, 'C');
utC = solution1.x(rxnEx) < 0 & Ematrix(:, eleC) > 0;
yieldC = solution1.f / abs(solution1.x(rxnEx(utC))' * Ematrix(utC, eleC)); 
fprintf('Symbiosis flux: %.4f g/C-mol\n', yieldC)

eleC = strcmp(elements, 'N');
utC = solution1.x(rxnEx) < 0 & Ematrix(:, eleC) > 0;
yieldN = solution1.f / abs(solution1.x(rxnEx(utC))' * Ematrix(utC, eleC)); 
fprintf('Symbiosis flux: %.4f g/N-mol\n', yieldN)

ratio = yieldC/yieldN;
fprintf('C/N: %.4f g/C-mol\n', ratio)

% C08
fprintf('Normalized symbiosis flux C08: %.4f\n', solution2.f);
rxnEx = find(sum(model.S ~= 0, 1) == 1);
[metEx, ~] = find(model.S(:, rxnEx));
[Masses, knownMasses, unknownElements, Ematrix, elements] = getMolecularMass(model.metFormulas(metEx), 0, true);
eleC = strcmp(elements, 'C');
utC = solution2.x(rxnEx) < 0 & Ematrix(:, eleC) > 0;
yieldC = solution2.f / abs(solution2.x(rxnEx(utC))' * Ematrix(utC, eleC)); 
fprintf('Symbiosis flux: %.4f g/C-mol\n', yieldC)

eleC = strcmp(elements, 'N');
utC = solution2.x(rxnEx) < 0 & Ematrix(:, eleC) > 0;
yieldN = solution2.f / abs(solution2.x(rxnEx(utC))' * Ematrix(utC, eleC)); 
fprintf('Symbiosis flux: %.4f g/N-mol\n', yieldN)

ratio = yieldC/yieldN;
fprintf('C/N: %.4f g/C-mol\n', ratio)

