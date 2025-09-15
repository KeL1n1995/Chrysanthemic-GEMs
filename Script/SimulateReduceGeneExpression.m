%% Data source: Fig.3 (D) (E)
%% The effects of simulated and experimental inhibition of ispA expression in CS5-1 strain on chrysanthemol titer.

clear
modelDBFile = 'Mains/Model/NewBL21DE3_CS5.mat';
load(modelDBFile);
% biomass
WT=5.038045655;% (mmol/g DW-hr)
CS5=1.706;% (mmol/g DW-hr)
CS_OH=0.0023;%(mmol/g DW-hr)

%% pre-set
List={'EX_colipa_e';
'EX_eca4colipa_e';
'EX_hacolipa_e';
'EX_o16a4colipa_e';
'EX_acolipa_e'};
model = changeRxnBounds(model, List,-1000*ones(5,1), 'l');
model = changeRxnBounds(model, List,0*ones(5,1), 'u');

model = changeRxnBounds(model, 'EX_glc__D_e',-20, 'l');
model = changeRxnBounds(model, 'EX_glc__D_e',0, 'u');

Model1=model;

%%
model = changeRxnBounds(model, 'DPMVD',-1000, 'l');
model = changeRxnBounds(model, 'DPMVD',1000, 'u');
model = changeRxnBounds(model, 'MEPCT',-1000, 'l');
model = changeRxnBounds(model, 'MEPCT',1000, 'u');

%% ADD ratio = MVA/MEP = 6461;
% Note:This value is calculated from transcript expression abundance in CS5 and WT
model = addRatioReaction(model, {'DPMVD' 'CDPMEK'}, [1 6461]); 

%% simulate

model = changeRxnBounds(model, 'EX_CS_OH',0.0023, 'l');
model = changeRxnBounds(model, 'EX_CS_OH',0.0023, 'u');

model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',CS5, 'l');
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',CS5, 'u');

model = changeObjective(model,'BIOMASS_Ec_iJO1366_WT_53p95M', 1);
FBA = optimizeCbModel(model,'max');
FBA.f;

g1=FBA.x((contains(model.rxns,'EX_glc__D_e')));%biomass
b1=FBA.x((contains(model.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M')));%biomass

af1=FBA.x((contains(model.rxns,'DMATT')));%grdp ispA
af2=FBA.x((contains(model.rxns,'GRTT')));%frdp ispA
af3=FBA.x((contains(model.rxns,'CSDP_syn'))); % TcCDS

af4=FBA.x((contains(model.rxns,'DPMVD')));%MVD1
af5=FBA.x((contains(model.rxns,'DMPPS')));%ispH
af6=FBA.x((contains(model.rxns,'IPDPS')));%ispH

input_f1=af4+af5+af6;

% The flux of the branching reaction involving DMPP cannot exceed the total flux of DMPP production 
% before ispA inhibition.
model = addCOBRAConstraints(model, {'DMATT','GRTT','CSDP_syn','OCTDPS','UDCPDPS','Lavandulol_dp_syn','DMAPP_DP'}, ...
input_f1*0.4, 'c', [1 1 1 1 1 1 1], 'dsense', 'L');


model = changeRxnBounds(model, 'EX_CS_OH',CS_OH, 'l');
model = changeRxnBounds(model, 'EX_CS_OH',1000, 'u');

%%
X=[0.5 1 1.5 2 2.5 3 3.5 4];% Note that when X = 0 here has no significance!!!
% The output before knockdown should be equal to the output in the initial CS5 strain. 
% Subsequent commands will reassign values. Z(1,1)=CS_OH;
inhibition=2.^X;
% The maximum inhibition fold of ISPA gene is 8
for i=1:length(X)
    x=inhibition(i);
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',CS5/x, 'l');% 
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',CS5, 'u');
model = changeRxnBounds(model, 'DMATT',af1/x, 'l');
model = changeRxnBounds(model, 'DMATT',af1/x, 'u');
model = changeRxnBounds(model, 'GRTT',af1/x, 'l');
model = changeRxnBounds(model, 'GRTT',af1/x, 'u');
% 
% model = changeObjective(model,'EX_CS_OH', 1);
% FBA = optimizeCbModel(model,'max');
% FBA.f;

model = changeObjective(model,{'BIOMASS_Ec_iJO1366_WT_53p95M','EX_CS_OH'}, 1);
FBA = optimizeCbModel(model,'max');
FBA.f;

Y(:,i)=FBA.x;
Z(1,i)=FBA.f;
end

%% Visualize the results
fluxBiomass=Y((contains(model.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M')),:);%biomass
fluxCS=Y((contains(model.rxns,'EX_CS_OH')),:);

fluxBiomass=[CS5 fluxBiomass];
fluxCS=[CS_OH fluxCS];

% Fig.3 (D) 
Ratio=fluxCS./fluxCS(1,1);
FC=[0 0.5 1 1.5 2 2.5 3 3.5 4]
plot(FC,Ratio);
xlabel(" -log2(FC(ispA's relative expression))");
ylabel("Chrysanthemol yield fold change");

