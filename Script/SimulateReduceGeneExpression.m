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

MEP.rxns={'DXPS';'DXPRIi';'MEPCT';'CDPMEK';'MECDPS';'MECDPDH2';'DMPPS';'IPDPS'};
MVA.rxns={'ACACT1r';'MHGS';'HMGCOAR';'MEVK1';'PMEVK';'DPMVD';'IPDDI';'DMATT';'GRTT'};

Model1=model;

%%
model = changeRxnBounds(model, 'DPMVD',-1000, 'l');
model = changeRxnBounds(model, 'DPMVD',1000, 'u');
model = changeRxnBounds(model, 'MEPCT',-1000, 'l');
model = changeRxnBounds(model, 'MEPCT',1000, 'u');

%% ADD ratio = MVA/MEP = 6461;
% Note:This value is calculated from transcript expression abundance in CS5 and WT
model = addRatioReaction(model, {'DPMVD' 'CDPMEK'}, [1 6461]); 

% % Set a minimum value for the target product
% model = changeRxnBounds(model, 'EX_CS_OH',0.0023, 'l');
% model = changeRxnBounds(model, 'EX_CS_OH',100, 'u');

model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',CS5, 'l');
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',100, 'u');

% model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',CS5, 'l');
% model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',CS5, 'u');

model = changeObjective(model,'EX_CS_OH', 1);
FBA = optimizeCbModel(model,'max');
max_CS=FBA.f;

N=FBA.x((contains(model.rxns,'CDPMEK')));%MEP  
M=FBA.x((contains(model.rxns,'DPMVD')));%MVA

model=Model1;

model = changeRxnBounds(model, 'CDPMEK',0, 'l');
model = changeRxnBounds(model, 'CDPMEK',N, 'u');

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
input_f1*1, 'c', [1 1 1 1 1 1 1], 'dsense', 'L');

%The biomass reaction is usually set to 1%-10% of maximum theoretical biomass yield 
%when running the following steps, to prevent solutions with not biomass formation:
%maximizing product formation

model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',b1*0.05, 'l');% origin set: 0.05%
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',b1, 'u');

model = changeRxnBounds(model, 'EX_CS_OH',CS_OH, 'l');
model = changeRxnBounds(model, 'EX_CS_OH',1000, 'u');


%%

% model = changeRxnBounds(model, 'OCTDPS',0, 'l');
% model = changeRxnBounds(model, 'OCTDPS',0, 'u');
% 
% model = changeRxnBounds(model, 'UDCPDPS',0, 'l');
% model = changeRxnBounds(model, 'UDCPDPS',0, 'u');

model = changeRxnBounds(model, 'Lavandulol_dp_syn',0, 'l');
model = changeRxnBounds(model, 'Lavandulol_dp_syn',0, 'u');

model = changeRxnBounds(model, 'DMAPP_DP',0, 'l');
model = changeRxnBounds(model, 'DMAPP_DP',0, 'u');


X=[0 0.5 1 1.5 2 2.5 3 3.5 4];% Note that when X = 0 here has no significance!!!
% The output before knockdown should be equal to the output in the initial CS5 strain. 
% Subsequent commands will reassign values. Z(1,1)=CS_OH;
inhibition=2.^X;
% The maximum inhibition fold of ISPA gene is 8
for i=1:length(X)
    x=inhibition(i);
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',b1/x, 'l');% 
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

Z(1,1)=CS_OH;


%% Visualize the results
fluxCS=Z;
fluxbiomass=Y((contains(model.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M')),:);%biomass
fluxbiomass(1,1)=CS5;
Production=fluxbiomass.*fluxCS;

Ratio=fluxCS./fluxCS(1,1);

fluxCS=Y((contains(model.rxns,'EX_CS_OH')),:);
fluxCS=[CS_OH fluxCS];

fluxBiomass=Y((contains(model.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M')),:);%biomass
fluxBiomass=[CS5 fluxBiomass];

Production=fluxBiomass.*fluxCS;

% Ratio=Production./Production(1,1);
Ratio=fluxCS./fluxCS(1,1);
allresult=[fluxCS;fluxBiomass;Production;Ratio];

plot(X,Ratio);

xlabel(" -log2(FC(ispA's relative expression))");
ylabel("chrysanthemol yield(mmol/g DW-hr)");


plot(Z)
R=Z./Z(1,1);
plot(R)
% 
% bar(X,Z);
% xlabel(" -log2(FC(ispA's relative expression))");
% ylabel("chrysanthemol yield(mmol/g DW-hr)");

% 
load('SelectRxns.mat')
idx=find(contains(model.rxns,SelectRxns));
SelectRxnsFluxs=Y(idx,:);
SelectRxnsNames=model.rxns(idx);











