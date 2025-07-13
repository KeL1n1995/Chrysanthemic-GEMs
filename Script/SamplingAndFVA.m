
WT=5.038;% (mmol/g DW-hr)
CS6=3.252;% (mmol/g DW-hr)
CS5=1.706;% (mmol/g DW-hr)
CS3=2.214;% (mmol/g DW-hr)
MEP.rxns={'DXPRIi';'MEPCT';'CDPMEK';'MECDPS';'MECDPDH2'};
MVA.rxns={'MHGS';'HMGCOAR';'MEVK1';'PMEVK';'DPMVD'};

%% simulate the WT strain(only MEP pathway)
model1=BL21_CS3;
% model1.c(:)=0;
model1.lb((contains(model1.rxns,MEP.rxns)))=0;
model1.ub((contains(model1.rxns,MEP.rxns)))=100;
model1.lb((contains(model1.rxns,MVA.rxns)))=0;
model1.ub((contains(model1.rxns,MVA.rxns)))=0;

model1 = changeRxnBounds(model1, 'EX_glc__D_e',-50, 'l');
model1 = changeRxnBounds(model1, 'EX_glc__D_e',0, 'u');

model1 = changeRxnBounds(model1, 'BIOMASS_Ec_iJO1366_WT_53p95M',WT, 'l');
model1 = changeRxnBounds(model1, 'BIOMASS_Ec_iJO1366_WT_53p95M',WT, 'u');

options.nStepsPerPoint = 200;
options.nPointsReturned = 10000;

[P_1, X1_1] =  sampleCbModel(model1, [], [], options);


%% simulate the CS3 strain(add MVA pathway)
% add ratio: MEP/MVA

% model2.c(:)=0;
model2 = addRatioReaction(BL21_CS3, {'DPMVD' 'CDPMEK'}, [1 5000]);

model2.lb((contains(model2.rxns,MEP.rxns)))=0;
model2.ub((contains(model2.rxns,MEP.rxns)))=0.1;
model2.lb((contains(model2.rxns,MVA.rxns)))=0;
model2.ub((contains(model2.rxns,MVA.rxns)))=100;

model2 = changeRxnBounds(model2, 'EX_glc__D_e',-50, 'l');
model2 = changeRxnBounds(model2, 'EX_glc__D_e',0, 'u');

model2 = changeRxnBounds(model2, 'BIOMASS_Ec_iJO1366_WT_53p95M',CS3, 'l');
model2 = changeRxnBounds(model2, 'BIOMASS_Ec_iJO1366_WT_53p95M',CS3, 'u');

options.nStepsPerPoint = 200;
options.nPointsReturned = 10000;
[P_2, X1_2] = sampleCbModel(model2, [], [], options);



nPts=200;
sample1=X1_1;
sample2=X1_2;
[sampleDiff, sampleRatio] = calcSampleDifference(sample1, sample2);

%% Calculate FVA
[minFlux1, maxFlux1] = fluxVariability(model1, 90);
[minFlux2, maxFlux2] = fluxVariability(model2, 90);

mets_bypass={'Geraniol_syn1' 'FPP_hydro' 'Farnesal_acid'};


minFlux1((contains(model1.rxns,mets_bypass)))
maxFlux1((contains(model1.rxns,mets_bypass)))

minFlux2((contains(model2.rxns,mets_bypass)))
maxFlux2((contains(model2.rxns,mets_bypass)))

List=[minFlux1, maxFlux1 minFlux2, maxFlux2];

% % Model = buildRxnEquations(modelOut2);
idx1 = ismember(model.rxns,'Geraniol_syn1'); % Geraniol
idx2 = ismember(model.rxns,'FPP_hydro'); % Farnesol
idx3 = ismember(model.rxns,'Farnesal_acid'); % Farnesal_acid

%% Visualization
clr_1 = [5,113,176]/255; %
clr_2 = [202,0,32]/255; %

subplot(1,3,1);
[y1, x1] = hist(X1_1(idx1, :), 20);
[y2, x2] = hist(X1_2(idx1, :), 20);
hold on;
plot(x1, y1, 'Color', clr_1, 'DisplayName', 'WT');
plot(x2, y2, 'Color', clr_2, 'DisplayName', 'CS3');
title('Geraniol');
xlabel('Flux (mmol/gCDW/h)');ylabel('Number of Samples');
legend('show');

subplot(1,3,2);
[y1, x1] = hist(X1_1(idx2, :), 20);
[y2, x2] = hist(X1_2(idx2, :), 20);
hold on;
plot(x1, y1, 'Color', clr_1, 'DisplayName', 'WT');
plot(x2, y2, 'Color', clr_2, 'DisplayName', 'CS3');
hold on;
title('Farnesol');
xlabel('Flux (mmol/gCDW/h)');ylabel('Number of Samples');
legend('show');

subplot(1,3,3);
[y1, x1] = hist(X1_1(idx3, :), 20);
[y2, x2] = hist(X1_2(idx3, :), 20);
hold on;
plot(x1, y1, 'Color', clr_1, 'DisplayName', 'WT');
plot(x2, y2, 'Color', clr_2, 'DisplayName', 'CS3');
title('Farnesal acid');
xlabel('Flux (mmol/gCDW/h)');ylabel('Number of Samples');
legend('show');


%% output data
data =X1_1 ; 
filename = 'SamplingData_WT.xlsx';
xlswrite(filename, data);
disp('Data has been written to Excel using xlswrite.');

data =X1_2 ; 
filename = 'SamplingData_CS3.xlsx';
xlswrite(filename, data);
disp('Data has been written to Excel using xlswrite.');

data=[X1_1(idx1, :);X1_1(idx2, :);X1_1(idx3, :);X1_2(idx1, :);X1_2(idx2, :);X1_2(idx3, :)];
filename = 'PartFlux.xlsx';
xlswrite(filename, data);
disp('Data has been written to Excel using xlswrite.');

bar(X1_1(idx1, :)); 
