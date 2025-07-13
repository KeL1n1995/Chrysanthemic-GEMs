
%% Data source: Fig.4ESimulated optimization of chrysanthemic acid production via enzyme flux manipulation.
clear
modelDBFile = 'Mains/Model/NewBL21DE3_CS5.mat';
load(modelDBFile);

%%  pre-set
List={'EX_colipa_e';
'EX_eca4colipa_e';
'EX_hacolipa_e';
'EX_o16a4colipa_e';
'EX_acolipa_e'};

model = changeRxnBounds(model, List,-1000*ones(5,1), 'l');
model = changeRxnBounds(model, List,0*ones(5,1), 'u');

model = changeRxnBounds(model, 'EX_glc__D_e',-20, 'l');
model = changeRxnBounds(model, 'EX_glc__D_e',0, 'u');

model = changeRxnBounds(model, 'DPMVD',-1000, 'l');
model = changeRxnBounds(model, 'DPMVD',1000, 'u');
model = changeRxnBounds(model, 'MEPCT',-1000, 'l');
model = changeRxnBounds(model, 'MEPCT',1000, 'u');

model = addRatioReaction(model, {'DPMVD' 'CDPMEK'}, [1 6000]); 
M_CS5=model;

%% Reconstriant BL21(E3)_CS6
model=M_CS5;

% CS6 chrysanthemum acid production, (AUQ44118.1) and ALDH1 (AUQ44119.1) 
% were introduced on the basis of CS5 and overexpressed.
% ADH2;Alcohol dehydrogenase 2
% (R,R)-chrysanthemol + NAD+ = (1R,3R)-chrysanthemal + H+ + NADH;   EC:1.2.1.5
model = addReaction(model,'ADH2','metaboliteList',...
    {'CSOH[c]','nad[c]','nadh[c]','CS=O[c]','h[c]'},...
    'stoichCoeffList',[-1 -1 1 1 1 ], 'reversible',false,'geneRule','TcADH2'); 
model = addReaction(model,'tran_CS=O','reactionFormula',...
    'CS=O[c] <=> CS=O[e]');                  
model = addReaction(model,'EX_CS=O','reactionFormula',...
    'CS=O[e] -> ');  

% ALDH1:Aldehyde dehydrogenase 1;
%(1R,3R)-chrysanthemal + H2O + NAD+ = (1R,3R)-chrysanthemate + 2 H+ + NADH;   EC:1.2.1.5
model = addReaction(model,'ALDH1','metaboliteList',...
    {'CS=O[c]','nad[c]','nadh[c]','CS-OOH[c]','h[c]'},...
    'stoichCoeffList',[-1 -1 1 1 1], 'reversible',false,'geneRule','TcALDH1'); 

model = addReaction(model,'tran_CS_OOH','reactionFormula',...
    'CS-OOH[c] <=> CS-OOH[e]');                  
model = addReaction(model,'EX_CS_OOH','reactionFormula',...
    'CS-OOH[e] -> ');  

% biomass
B_CS6=1.238940474;% (mmol/g DW-hr)

model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',B_CS6, 'l');
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',B_CS6, 'u');

M_CS6=model;

%%
model=M_CS6;


model = addRatioReaction(model, {'ADH2' 'ALDH1'}, [1 6.65]); % flux_ADH2:flux_ALDH1=6.650052833
model = changeObjective(model,'EX_CS_OOH', 1);
FBA = optimizeCbModel(model,'max');
FBA.f;
max_CS6=FBA.f;

f1=FBA.x((contains(model.rxns,'ADH2')));%% ADH2;Alcohol dehydrogenase 2
f2=FBA.x((contains(model.rxns,'ALDH1')));%% ALDH1:Aldehyde dehydrogenase 1;

%% Simulate


X=1:1:5;
inhibition=2.^X;

n=0;
flux1=[];
model=M_CS6;
for i=1:1:5
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',B_CS6*0.2, 'l');
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',B_CS6, 'u');
model = changeRxnBounds(model, 'ADH2',0, 'l');
model = changeRxnBounds(model, 'ADH2',f1, 'u');    
model = changeRxnBounds(model, 'ALDH1',0, 'l');
model = changeRxnBounds(model, 'ALDH1',f2*inhibition(i), 'u');
model = changeObjective(model,'EX_CS_OOH', 1);
FBA = optimizeCbModel(model,'max');
FBA.f;
n=n+1;
flux1(n,1)=FBA.f;
end


n=0;
flux2=[];
model=M_CS6;
for i=1:1:5
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',B_CS6*0.2, 'l');
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',B_CS6, 'u');
model = changeRxnBounds(model, 'ADH2',0, 'l');
model = changeRxnBounds(model, 'ADH2',f1*inhibition(i), 'u');
model = changeRxnBounds(model, 'ALDH1',0, 'l');
model = changeRxnBounds(model, 'ALDH1',f2, 'u');
model = changeObjective(model,'EX_CS_OOH', 1);
FBA = optimizeCbModel(model,'max');
FBA.f;
n=n+1;
flux2(n,1)=FBA.f;
end

n=0;
flux3=[];
model=M_CS6;
for i=1:1:5
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',B_CS6*0.2, 'l');
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',B_CS6, 'u');
model = changeRxnBounds(model, 'ADH2',0, 'l');
model = changeRxnBounds(model, 'ADH2',f1*inhibition(i), 'u');
model = changeRxnBounds(model, 'ALDH1',0, 'l');
model = changeRxnBounds(model, 'ALDH1',f2*inhibition(i), 'u');
model = changeObjective(model,'EX_CS_OOH', 1);
FBA = optimizeCbModel(model,'max');
FBA.f;
n=n+1;
flux3(n,1)=FBA.f;
end

CombineResults=[flux1 flux2 flux3];


%% Optional: --- Visualization ---

% The x-axis values are now the geometric series
x_values = 1:1:5 ;

figure; % Create a new figure window for the plot

% Plot flux1 (varying ALDH1 only)
plot(x_values, flux1, 'b-o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'DisplayName', 'Flux when only ALDH1 varies');
hold on; % Keep the current plot active to add more lines

% Plot flux2 (varying ADH2 only)
plot(x_values, flux2, 'r-s', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'DisplayName', 'Flux when only ADH2 varies');

% Plot flux3 (varying both ADH2 and ALDH1)
plot(x_values, flux3, 'g-^', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'DisplayName', 'Flux when both ADH2 & ALDH1 vary');

hold off; % Release the plot

% Customize the Plot
title('Max Chrysanthemic Yield Flux under Different Constraint Scenarios');
xlabel(" log2(Multiplier for Base Flux)");
ylabel('Maximum Chrysanthemic Flux (mmol/gDW/h)');
legend('show', 'Location', 'best');
grid on;
