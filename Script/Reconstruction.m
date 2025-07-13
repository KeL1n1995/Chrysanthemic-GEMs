
clear

modelDBFile = 'Mains/Model/NewBL21DE3_WT.mat';
load(modelDBFile);


% model = readCbModel('NewBL21DE3_WT','fileType','SBML');

[modelout, ~, ~] =checkDuplicateRxn(model,'FR');

model=modelout;

str = model.mets; % Original text
expression = '_([a-z])$'; % regular expression that matches underscores and a lowercase letter at the end
replace = '[$1]'; % alts text, enclosing matching letters in square brackets
newStr = regexprep(str, expression, replace); % uses the regexprep function
model.mets=newStr;

[missingMets, presentMets, coupledMets, missingCofs, presentCofs] = ...
    biomassPrecursorCheck(model,true,true);


%% Gap filling bypass reaction
 %% lavandulol
% EC Number: 2.5.1.69
% 2 prenyl diphosphate(DMPP)  lavandulyl diphosphate + diphosphate(EC:2.5.1.69 FDS-5:lavandulyl diphosphate synthase)
model = addReaction(model,'Lavandulol_dp_syn','metaboliteList',{'dmpp[c]','ppi[c]','lavandulol-dp[c]'},'stoichCoeffList',[-2 1 1], 'reversible',false,'geneRule', 'FDS-5');
model= ReviseMetsFormula(model,'lavandulol-dp[c]', {'C10H20O7P2'});% Lavandulyl diphosphate

% lavandulyl diphosphate + H2O  lavandulol + diphosphate [no EC number assigned]
model = addReaction(model,'Lavandulol_syn','metaboliteList',...
    {'lavandulol-dp[c]','h2o[c]','lavandulol[c]','ppi[c]'},'stoichCoeffList',[-1 -1 1 1], 'reversible',false,'geneRule', '');
model = addReaction(model,'Trans_lavandulol','reactionFormula',...
    'lavandulol[c] <=> lavandulol[e]');                  
model = addReaction(model,'EX_lavandulol','reactionFormula',...
    'lavandulol[e] -> ');  

model= ReviseMetsFormula(model,'lavandulol[c]', {'C10H18O'});% 
model= ReviseMetsFormula(model,'lavandulol[e]', {'C10H18O'});% 

%% Geraniol
% Geranyl diphosphate + H2O <=> Geraniol + Diphosphate EC:3.1.7.3///3.1.7.11- geranyl diphosphate diphosphatase

model = addReaction(model,'Geraniol_syn1','metaboliteList',...
    {'grdp[c]','h2o[c]','Geraniol[c]','ppi[c]'},'stoichCoeffList',[-1 -1 1 1], 'reversible',false,'geneRule','GES');
% % Geranyl phosphate + H2O <=> Geraniol + pi[c]  EC:3.1.7.3///3.1.7.11
% model = addReaction(model,'Geraniol_syn2','metaboliteList',...
%     {'grmp[c]','h2o[c]','Geraniol[c]','pi[c]'},'stoichCoeffList',[-1 -1 1 1], 'reversible',false,'geneRule', '');
model = addReaction(model,'Trans_Geraniol','reactionFormula',...
    'Geraniol[c] <=> Geraniol[e]');                  
model = addReaction(model,'EX_Geraniol','reactionFormula',...
    'Geraniol[e] -> ');  

model= ReviseMetsFormula(model,'Geraniol[c]', {'C10H18O'});% 
model= ReviseMetsFormula(model,'Geraniol[e]', {'C10H18O'});% 

%% Farnesol
%Phosphatases hydrolyze FPP to farnesol.[gene=ybjG] [locus_tag=K5T46_RS04090][gene=pgpB] [locus_tag=K5T46_RS11165]	
model = addReaction(model,'FPP_hydro','metaboliteList',...
    {'frdp[c]','h2o[c]','Farnesol[c]','ppi[c]'},'stoichCoeffList',[-1 -1 1 1], 'reversible',false,'geneRule','K5T46_RS04090 and K5T46_RS11165');

model = addReaction(model,'Trans_Farnesol','reactionFormula',...
    'Farnesol[c] <=> Farnesol[e] ');                  
model = addReaction(model,'EX_Farnesol','reactionFormula',...
    'Farnesol[e] -> ');  

model= ReviseMetsFormula(model,'Farnesol[c]', {'C15H26O'});% 
model= ReviseMetsFormula(model,'Farnesol[e]', {'C15H26O'});% 


% NAME: (2E,6E)-farnesol:NAD+ 1-oxidoreductase
% [protein=NAD-dependent epimerase/dehydratase family protein]
% 	ENZYME: 1.2.1.3

model = addReaction(model,'Farnesol_oxidase','metaboliteList',...
    {'Farnesol[c]','nad[c]','Farnesal[c]','nadh[c]','h[c]'},'stoichCoeffList',[-1 -1 1 1 1], 'reversible',false,'geneRule','K5T46_RS00140');
model = addReaction(model,'Trans_Farnesal','reactionFormula',...
    'Farnesal[c] <=> Farnesal[e] ');                  
model = addReaction(model,'EX_Farnesal','reactionFormula',...
    'Farnesal[e] -> ');  

model= ReviseMetsFormula(model,'Farnesal[c]', {'C15H24O'});% 
model= ReviseMetsFormula(model,'Farnesal[e]', {'C15H24O'});% 

%farnesal:NAD+ oxidoreductase
% R08146
% [gene=astD] [locus_tag=K5T46_RS02595] 
% [locus_tag=K5T46_RS00140] [protein=NAD-dependent epimerase/dehydratase family protein] 
model = addReaction(model,'Farnesal_acid','metaboliteList',...
    {'Farnesal[c]','nad[c]','h2o[c]','Farnesoic_acid[c]','h[c]','nadh[c]'},...
    'stoichCoeffList',[-1 -1 -1 1 1 1], 'reversible',true,'geneRule','K5T46_RS00140 and K5T46_RS02595');

model = addReaction(model,'Trans_Farnesoic_acid','reactionFormula',...
    'Farnesoic_acid[c] <=> Farnesoic_acid[e] ');                  
model = addReaction(model,'EX_Farnesoic_acid','reactionFormula',...
    'Farnesoic_acid[e] -> ');  

model= ReviseMetsFormula(model,'Farnesoic_acid[c]', {'C15H24O2'});% 
model= ReviseMetsFormula(model,'Farnesoic_acid[e]', {'C15H24O2'});% 

%%  AddBypassPathwayReactions
        %%%%%% Geranyl diphosphate dephosphorylation NudI      
model = addReaction(model,'NudI','metaboliteList',...
    {'grdp[c]','h2o[c]','geranyl-P[c]','pi[c]','h[c]'},...
    'stoichCoeffList',[-1 -1 1 1 1], 'reversible',false,'geneRule','K5T46_RS13295');%                 
model = addReaction(model,'geranyl-P_trans','reactionFormula',...
    'geranyl-P[c] <=> geranyl-P[e] ');                  
model = addReaction(model,'EX_geranyl-P','reactionFormula',...
    'geranyl-P[e] -> ');  

model= ReviseMetsFormula(model,'geranyl-P[c]', {'C10H17O4P'});% geranyl monophosphate
model= ReviseMetsFormula(model,'geranyl-P[e]', {'C10H17O4P'});% geranyl monophosphate


        %%%%%% DMAPP dephosphorylation NudJ or NudF                
model = addReaction(model,'DMAPP_DP','metaboliteList',...
    {'dmpp[c]','h2o[c]','dmap[c]','pi[c]','h[c]'},...
    'stoichCoeffList',[-1 -1 1 1 1], 'reversible',false,'geneRule','K5T46_RS00410 or K5T46_RS12370');                
                
model = addReaction(model,'dmap_trans','reactionFormula',...
    'dmap[c] <=> dmap[e] ');                 
model = addReaction(model,'EX_dmap','reactionFormula',...
    'dmap[e] -> '); 

        %%%%%% prenyl diphosphate (DMPP)+ H2O = prenol + diphosphate / EC 3.1.7.1 / R03541

model= ReviseMetsFormula(model,'dmap[c]', {'C5H11O4P'});% Dimethylallyl phosphate
model= ReviseMetsFormula(model,'dmap[e]', {'C5H11O4P'});% Dimethylallyl phosphate

model = addReaction(model,'YAJO','metaboliteList',...
    {'AH2[c]','ru5p__D[c]','dxyl5p[c]','A[c]','h2o[c]'},...
    'stoichCoeffList',[-1 -1 1 1 1], 'reversible',false,'geneRule','K5T46_RS15850');                               
model = addReaction(model,'AH2_trans','reactionFormula',...
    'AH2[c] <=> AH2[e] ');                  
model = addReaction(model,'EX_AH2','reactionFormula',...
    'AH2[e] -> '); 
model = addReaction(model,'A_trans','reactionFormula',...
    'A[c] <=> A[e] ');                 
model = addReaction(model,'EX_A','reactionFormula',...
    'A[e] -> ');            

model= ReviseMetsFormula(model,'AH2[c]', {'RH2'});%  
model= ReviseMetsFormula(model,'AH2[e]', {'RH2'});%  
model= ReviseMetsFormula(model,'A[c]', {'R'});%  
model= ReviseMetsFormula(model,'A[e]', {'R'});% 

model = updateGenes(model);  
model = buildRxnGeneMat(model);    
model.proteins=model.genes;
model.geneNames=model.genes;    
BL21_WT=model;

%% ADD MVA PATHWAY
% CS3 terpene chassis, introduced and overexpressed atoB, tHMGR, HMGS, ERG12, ERG8, MVD1 idi genes.
% Acetyl-CoA acetyltransferase(atoB)% EC:2.3.1.9 
% 2 accoa[c]  <=> coa[c] + aacoa[c] %add  gene
model = changeGeneAssociation(model,'ACACT1r','K5T46_RS13160 or K5T46_RS12775 or K5T46_RS09705 or TcatoB');

% Hydroxymethylglutaryl-CoA synthase(ERG13)% EC:2.3.3.10
% acetoacetyl-CoA + acetyl-CoA + H2O = (3S)-hydroxy-3-methylglutaryl-CoA + CoA + H+
model = addReaction(model,'HMGS','metaboliteList',...
    {'aacoa[c]','accoa[c]','h2o[c]','h[c]','coa[c]','hmgcoa[c]'},...
    'stoichCoeffList',[-1 -1 -1 1 1 1], 'reversible',false,'geneRule','TcHMGS');      

%tHMGR %
% HMGCOAR: Hydroxymethylglutaryl CoA reductase  (tHMGR)
% (R)-mevalonate + CoA + 2 NAD(+) <=> (3S)-hydroxy-3-methylglutaryl-CoA + 2 H(+) + 2 NADH(EC 1.1.1.88)
% (R)-mevalonate + CoA + 2 NAD(+) <=> (3S)-hydroxy-3-methylglutaryl-CoA + 2 H(+) + 2 NADPH(EC 1.1.1.34)

model = addReaction(model,'HMGCOAR','metaboliteList',...
    {'h[c]','nadph[c]','hmgcoa[c]','coa[c]','nadp[c]','mev__R[c]'},...
    'stoichCoeffList',[-2 -2 -1 1 2 1], 'reversible',false,'geneRule','TctHMGR');

% MEVK1x: Mevalonate kinase atp
model = addReaction(model,'MEVK1','metaboliteList',...
    {'atp[c]','mev__R[c]','adp[c]','h[c]','5pmev[c]'},...
    'stoichCoeffList',[-1 -1 1 1 1], 'reversible',false,'geneRule','TcERG12');    

% PMEVK: Phosphomevalonate kinase
model = addReaction(model,'PMEVK','metaboliteList',...
    {'atp[c]','5pmev[c]','adp[c]','5dpmev[c]'},...
    'stoichCoeffList',[-1 -1 1 1], 'reversible',false,'geneRule','TcERG8'); 

% DPMVD: Diphosphomevalonate decarboxylase
model = addReaction(model,'DPMVD','metaboliteList',...
    {'atp[c]','5dpmev[c]','adp[c]','co2[c]','ipdp[c]','pi[c]'},...
    'stoichCoeffList',[-1 -1 1 1 1 1], 'reversible',false,'geneRule','TcMVD1'); 
% idi {'IPDDI'} (Already exists in the model)

model = changeGeneAssociation(model,'IPDDI','K5T46_RS05990 or Tcidi');

model= ReviseMetsFormula(model,'hmgcoa[c]', {'C27H39N7O20P3S'});% 
model= ReviseMetsFormula(model,'mev__R[c]', {'C6H11O4'});% 
model= ReviseMetsFormula(model,'5pmev[c]', {'C6H10O7P'});% 
model= ReviseMetsFormula(model,'5dpmev[c]', {'C6H11O10P2'});% 

BL21_CS3=model;

%% Reconstriant BL21(E3)_CS5 
% CS4 production of inulin, import and overexpress CDS gene (EC/2.5.1.67) on the basis of CS3.
% 2 dimethylallyl diphosphate <=> (R,R)-chrysanthemyl diphosphate(CDP) + diphosphate
% TcCDS : Chrysanthemyl diphosphate synthase 
% 
model = addReaction(model,'CSDP_syn','metaboliteList',...
    {'dmpp[c]','ppi[c]','CSDP[c]'},...
    'stoichCoeffList',[-2 1 1], 'reversible',false,'geneRule','TcCDS'); 

%% CS5 methrin production, introduced and overexpressed Nudix1 gene (QIJ31369.1) on the basis of CS4.
% (R,R)-chrysanthemyl diphosphate(-3) + h2o = (R,R)-chrysanthemol + diphosphate
model = addReaction(model,'NUDIX1','metaboliteList',...
    {'CSDP[c]','h2o[c]','ppi[c]','CSOH[c]'},...
    'stoichCoeffList',[-1 -1 1 1], 'reversible',false,'geneRule','TcNudix1'); 

model = addReaction(model,'tran_CS_OH','reactionFormula',...
    'CSOH[c] <=> CSOH[e]');                  
model = addReaction(model,'EX_CS_OH','reactionFormula',...
    'CSOH[e] -> ');


model= ReviseMetsFormula(model,'CSDP[c]', {'C10H20O7P2'});%  Chrysanthemyl diphosphate
model= ReviseMetsFormula(model,'CSOH[c]', {'C10H18O'});%  (R,R)-chrysanthemol
model= ReviseMetsFormula(model,'CSOH[e]', {'C10H18O'});%  (R,R)-chrysanthemol


model = updateGenes(model);  
model = buildRxnGeneMat(model);    
model.proteins=model.genes;
model.geneNames=model.genes;

[modelOut, removedRxnInd, keptRxnInd]=checkDuplicateRxn(model,'FR');

BL21_CS5=modelOut;
model=BL21_CS5;
m = length(model.subSystems); %
cellArray = cell(m, 1); % 
cellArray(:) = {''}; % 
model.subSystems=cellArray;
model.modelID='NewBL21DE3_CS5';

%% Reconstriant BL21(E3)_CS6
% CS6 chrysanthemum acid production, (AUQ44118.1) and ALDH1 (AUQ44119.1) 
% were introduced on the basis of CS5 and overexpressed.
% ADH2;Alcohol dehydrogenase 2
% (R,R)-chrysanthemol + NAD+ = (1R,3R)-chrysanthemal + H+ + NADH;   EC:1.2.1.5
model = addReaction(model,'ADH2','metaboliteList',...
    {'CSOH[c]','nad[c]','nadh[c]','CS=O[c]','h[c]'},...
    'stoichCoeffList',[-1 -1 1 1 1 ], 'reversible',false,'geneRule','TcADH2'); 

%ADH1:Aldehyde dehydrogenase 1;
%(1R,3R)-chrysanthemal + H2O + NAD+ = (1R,3R)-chrysanthemate + 2 H+ + NADH;   EC:1.2.1.5
model = addReaction(model,'ALDH1','metaboliteList',...
    {'CS=O[c]','nad[c]','nadh[c]','CS-OOH[c]','h[c]'},...
    'stoichCoeffList',[-1 -1 1 1 2], 'reversible',false,'geneRule','TcALDH1'); 

model = addReaction(model,'tran_CS_OOH','reactionFormula',...
    'CS-OOH[c] <=> CS-OOH[e]');                  
model = addReaction(model,'EX_CS_OOH','reactionFormula',...
    'CS-OOH[e] -> ');  










