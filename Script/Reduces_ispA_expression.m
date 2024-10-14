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

model = changeRxnBounds(model, 'EX_glc__D_e',-50, 'l');
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

% Set a minimum value for the target product
model = changeRxnBounds(model, 'EX_CS_OH',0.0023, 'l');
model = changeRxnBounds(model, 'EX_CS_OH',100, 'u');

model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',CS5, 'l');
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',100, 'u');

model = changeObjective(model,'EX_CS_OH', 1);
FBA = optimizeCbModel(model,'max');
FBA.f;

N=FBA.x((contains(model.rxns,'CDPMEK')));%MEP  
M=FBA.x((contains(model.rxns,'DPMVD')));%MVA

model=Model1;
model = changeRxnBounds(model, 'CDPMEK',0, 'l');
model = changeRxnBounds(model, 'CDPMEK',N, 'u');
%%

model = changeRxnBounds(model, 'EX_CS_OH',0.0023, 'l');
model = changeRxnBounds(model, 'EX_CS_OH',0.0023, 'u');

model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',0, 'l');
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',10, 'u');

model = changeObjective(model,'BIOMASS_Ec_iJO1366_WT_53p95M', 1);
FBA = optimizeCbModel(model,'max');
FBA.f;

g1=FBA.x((contains(model.rxns,'EX_glc__D_e')));%biomass
b1=FBA.x((contains(model.rxns,'BIOMASS_Ec_iJO1366_WT_53p95M')));%biomass
f1=FBA.x((contains(model.rxns,'DMATT')));%grdp
f2=FBA.x((contains(model.rxns,'GRTT')));%frdp
f3=FBA.x((contains(model.rxns,'DPMVD')));
f4=FBA.x((contains(model.rxns,'CSDP_syn')));

%The biomass reaction is usually set to 1%-10% of maximum theoretical biomass yield 
%when running the following steps, to prevent solutions with not biomass formation:
%maximizing product formation
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',0.15*b1, 'l');
model = changeRxnBounds(model, 'BIOMASS_Ec_iJO1366_WT_53p95M',b1, 'u');
model = changeRxnBounds(model, 'EX_CS_OH',0.0023, 'l');
model = changeRxnBounds(model, 'EX_CS_OH',30, 'u');

% The flux of the branching reaction involving DMPP cannot exceed the total flux of DMPP production 
% before ispA inhibition.
model = addCOBRAConstraints(model, {'DMATT','DPMVD'}, f3+f2+f1, 'c', [1 1], 'dsense', 'L');

% The maximum inhibition fold of ISPA gene is 5

X=[0 0.5 1 1.5 2 2.5];
inhibition=2.^X;
for i=1:6
    x=inhibition(i);
model = changeRxnBounds(model, 'DMATT',f1/x, 'l');
model = changeRxnBounds(model, 'DMATT',f1/x, 'u');
model = changeRxnBounds(model, 'GRTT',f1/x, 'l');
model = changeRxnBounds(model, 'GRTT',f1/x, 'u');
model = changeObjective(model,'EX_CS_OH', 1);
FBA = optimizeCbModel(model,'max');
FBA.f;
Y(:,i)=FBA.x;
Z(1,i)=FBA.f;
end

scatter(X,Z,'fill');

bar(X,Z);

xlabel(" -log2(FC(ispA's relative expression))");
ylabel("chrysanthemol yield(mmol/g DW-hr)");


