 
%% balance Chemical elements 
% pre-set 
modelDBFile = 'BiGG_universal.mat';
load(modelDBFile);
modelDB = model;

load('Database/Compounds.mat')
MetInfor=table2cell(Compounds);

load('Model/NewBL21DE3.mat');

model.rev = zeros(length(model.rxns),1);
model.rev(intersect(find(model.lb < 0), find(model.ub > 0))) = 1;
model = changeObjective (model, 'BIOMASS_Ec_iJO1366_WT_53p95M');

model= removeRxns(model, 'BIOMASS_Ec_iJO1366_core_53p95M');% remove the redundancy biomass rxns
model= removeRxns(model, 'LIPAcore');% remove the redundancy  rxns
model= removeRxns(model, 'EDTXS_core');% remove the redundancy  rxns
model= removeRxns(model, 'EDTXSCOF');% remove the redundancy  rxns


model = addMetabolite(model,'colipaA[c]','metFormula','C161H260N2O88P4','Charge',0);
model = addReaction(model,'GALR1TRA2','gagggicolipaAR1[c] + udpgal[c] -> h[c] + udp[c] + colipaA[c]');
model = addReaction(model,'O6GLCT1','udpg[c] + gammagund[c] -> udp[c] + o6aund[c] + h[c] + nh4[c]');
model = addReaction(model,'RIBCTR','reactionFormula','ppi[c] + cdprbt[c] <=> ctp[c] + h[c] + rbt5p[c]');
model = addReaction(model,'ACGAL6PISO','reactionFormula','h2o[c] + acgal6p[c] -> ac[c] + galam6p[c] + h[c]');


model.metFormulas = replace(model.metFormulas,'P0','P');% Replace 'P0' in phosphorus-containing metabolite formulas

rxnsToIgnore = ones(length(model.rxns),1); % Do not consider exchange, source, sink, and biomass reactions since they are inherently imbalanced
rxnsToIgnore((findExcRxns(model))) = 0;
rxnsToIgnore((model.c==1)) = 0;
model.SIntRxnBool = logical(rxnsToIgnore);

Met_Index01=find( contains(model.mets,'LptA[p]'));
cellArray=model.metFormulas(Met_Index01);
cellArrayWithR = strcat(cellArray, 'R');
model.metFormulas(Met_Index01)=cellArrayWithR;

[~,Met_Index1] = ismember('colipa_LptA[p]',model.mets);
[~,Met_Index2] = ismember('colipa[p]',model.mets);
model.metFormulas(Met_Index1)=strcat(model.metFormulas(Met_Index2), 'R');

model= ReviseMetsFormula(model,'o6a4colipa[p]', {'C312H523N6O200P4'});% o16a4colipa
model= ReviseMetsFormula(model,'o6a4colipa[e]', {'C312H523N6O200P4'});% o16a4colipa
model= ReviseMetsFormula(model,'o6a4colipa_LptA[p]', {'C312H523N6O200P4R'});% o16a4colipa

model= ReviseMetsFormula(model,'o6a4und[p]', {'C191H310N4O107P2'});% o16a4und 
model= ReviseMetsFormula(model,'o6aund[p]',{'C89H145NO32P2'});% o16aund 
model= ReviseMetsFormula(model,'o6aund[c]',{'C89H145NO32P2'});% o16aund 
model= ReviseMetsFormula(model,'o6a2und[p]', {'C123H200N2O57P2'});% o16a2und 
model= ReviseMetsFormula(model,'o6a3und[p]', {'C157H255N3O82P2'});% o16a3und 
model= ReviseMetsFormula(model,'icolipaBR1[c]', {'C129H206N2O63P4'});% 
model= ReviseMetsFormula(model,'hphhlipaB[c]', {'C129H207N2O60P3'});% 
model= ReviseMetsFormula(model,'phhlipaB[c]', {'C122H195N2O54P3'});% 
model= ReviseMetsFormula(model,'colipaF[p]', {'C176H302N2O103P5'});%  
model= ReviseMetsFormula(model,'lipidA_core[p]', {'C68H126N2O23P2'});%  
model= ReviseMetsFormula(model,'hhlipaA[c]', {'C124H220N2O51P2'});% 
model= ReviseMetsFormula(model,'hhlipaB[c]', {'C122H216N2O51P2'});% 
model= ReviseMetsFormula(model,'phhlipaA[c]', {'C124H219N2O54P3  '});% 
model= ReviseMetsFormula(model,'phhlipaB[c]', {'C122H215N2O54P3 '});% 
model= ReviseMetsFormula(model,'hphhlipaB[c]', {'C129H227N2O60P3'});% 
model= ReviseMetsFormula(model,'icolipaAR1[c]', {'C131H230N2O63P4'});% 
model= ReviseMetsFormula(model,'icolipaBR1[c]', {'C129H226N2O63P4'});% 
model= ReviseMetsFormula(model,'gicolipaAR1[c]', {'C137H240N2O68P4'});% 
model= ReviseMetsFormula(model,'gicolipaBR1[c]', {'C135H236N2O68P4'});% 
model= ReviseMetsFormula(model,'ggicolipaAR1[c]', {'C143H250N2O73P4'});% 
model= ReviseMetsFormula(model,'ggicolipaBR1[c]', {'C141H246N2O73P4'});% 
model= ReviseMetsFormula(model,'gggicolipaAR1[c]', {'C149H260N2O78P4'});% 
model= ReviseMetsFormula(model,'gggicolipaBR1[c]', {'C147H256N2O78P4'});% 
model= ReviseMetsFormula(model,'gagggicolipaAR1[c]', {'C155H270N2O83P4'});% 
model= ReviseMetsFormula(model,'gagggicolipaBR1[c]', {'C153H266N2O83P4'});% 
model= ReviseMetsFormula(model,'colipaB[c]', {'C159H276N2O88P4'});% 
model= ReviseMetsFormula(model,'colipaA[c]', {'C161H280N2O88P4'});% 
model= ReviseMetsFormula(model,'colipaB[p]', {'C159H276N2O88P4'});% 

[~,imBalancedMass,~,~,~,~,~] = checkMassChargeBalance(model,-1);
imbalancedRxnsMass = setdiff(find(~cellfun(@isempty,imBalancedMass)),find(rxnsToIgnore == 0));

rxnList_im01=model.rxns(imbalancedRxnsMass);

Fmodel=model;
Fmodel.mets=model.metFormulas;
Formulas_im1=printRxnFormula(model,'rxnAbbrList',rxnList_im01);
Formulas_im2=printRxnFormula(Fmodel,'rxnAbbrList',rxnList_im01);
imMass=imBalancedMass(imbalancedRxnsMass);

testR01=[Formulas_im1 Formulas_im2 imMass];


for r = 1:length(imbalancedRxnsMass)
        imbalancedElementsCurr = imBalancedMass{imbalancedRxnsMass(r)};
        s = split(imbalancedElementsCurr);
%         Str{r,1}= s;
        if length(s) == 2 && strcmp(s{2},'H')
            coeff = str2double(s{1});
            if  abs(coeff) <=5
            model.S(find(ismember(model.mets,'h[c]')),imbalancedRxnsMass(r)) = model.S(find(ismember(model.mets,'h[c]')),imbalancedRxnsMass(r)) - coeff;
            end
        end
end


strings=model.mets;
modifiedStrings = regexprep(strings, '\[([a-zA-Z])\]$', '_$1');
model.metSBOTerms = repmat({'SBO:0000247'}, length(model.mets), 1);  

    if ~isfield(model,'metCompSymbol')
        metCompSymbol = cell(length(model.mets),1);
        for m = 1:length(modifiedStrings)
            met = model.mets{m};
            symbol = met(end);
            metCompSymbol{m} = symbol;
        end
        model.metCompSymbol = metCompSymbol;
    end

% define compartment
model.comps = {'c';'e';'p'};  
% Name the compartment
model.compNames = {'cytosol'; 'extracellular';'periplasm'};
% Ensures each metabolite is assigned to the correct compartment
for i = 1:length(model.mets)
    if endsWith(model.mets{i}, '[c]')
        model.metComps(i,1) = find(strcmp(model.comps, 'c'));
    elseif endsWith(model.mets{i}, '[e]')
        model.metComps(i,1) = find(strcmp(model.comps, 'e'));
    elseif endsWith(model.mets{i}, '[p]')
        model.metComps(i,1) = find(strcmp(model.comps, 'p'));
    else
        warning('Metabolite %s does not match any known compartments.', model.mets{i});
    end
end    
    
% model = addSBOterms(model);% add rxns SBOTerms  annotateSBOTerms(model)

m = length(model.subSystems); %
cellArray = cell(m, 1); % 
cellArray(:) = {''}; % 
model.subSystems=cellArray;

model.modelID='NewBL21DE3_WT';

outmodel = writeCbModel(model, 'format','mat', 'fileName', 'Model/NewBL21DE3_WT.mat');

%% Optional : get memote report
% sbmlModel = writeSBML (model, 'NewBL21DE3_WT');
% command = 'memote report snapshot --filename "Memote_BL21DE3_WT.html" NewBL21DE3_WT.xml'; % 
% status = system(command);




