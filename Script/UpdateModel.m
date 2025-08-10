%% read genome data & input model file

[id1,~] = fastaread('Database/CDS_CP_001509.3.fna');% gene id format:ECD_xxxxx;Old model gene name 
[id2,~]= fastaread('Database/CDS_NC_012971.2.fna');% gene id format:ECD_RSxxxx;Old model gene name 
[id3,~]= fastaread('Database/CDS_NZ_CP081489.1.fna'); % gene id format:K5T46_RSxxxx; New model gene name 

ori_model = readSBML('Model/iEC1356_Bl21DE3.xml',1000);% input orgin model
ori_model.rxns=replace(ori_model.rxns,'_copy2','');
ori_model.rxns=replace(ori_model.rxns,'_copy1','');
[ori_model, ~, ~, ~] = removeGenesFromModel(ori_model,{'s0001' 'deleted'});

[modelOut, removedRxnInd, keptRxnInd] = checkDuplicateRxn(ori_model, 'FR');

model=modelOut;
model.rev = zeros(length(model.rxns),1);
model.rev(intersect(find(model.lb < 0), find(model.ub > 0))) = 1;

Origin_gene=model.genes;

%% Extract and classify gene ID information
info1=ExtractProteinInformation(id1);
info2=ExtractProteinInformation(id2);
info3=ExtractProteinInformation(id3);

%% map gene id
% input BLASTp result
table01=readcell('Database/BLASTp_CP2NZ.txt');
table02=readcell('Database/BLASTp_CP2NC.txt');
table03=readcell('Database/BLASTp_NC2NZ.txt');

% table01 = readtable('Database/BLASTp_CP2NZ.txt', 'FileType', 'text');% CP map NZ
% table02 = readtable('Database/BLASTp_CP2NC.txt', 'FileType', 'text');% CP map NC
% table03 = readtable('Database/BLASTp_NC2NZ.txt', 'FileType', 'text');% NC map NZ


S1_ID = info1.ID;
S2_ID = info2.ID;
S3_ID = info3.ID;

% CP map NZ
A1 = table01(:, 1);
A2 = table01(:, 2);
MapList01 = cell(length(A1), 2);

[~, Lia1] = ismember(A1, S1_ID);
[~, Lia2] = ismember(A2, S3_ID);

MapList01(:, 1) = info1.locus_tag(Lia1);
MapList01(:, 2) = info1.protein(Lia1);
MapList01(:, 3) = info1.gene(Lia1);
MapList01(:, 4) = info3.locus_tag(Lia2);
MapList01(:, 5) = info3.protein(Lia2);
MapList01(:, 6) = info3.gene(Lia2);
MapList01(:, 7)=table01(:,11);%E value

% NC map NZ
C1 = table03(:, 1);
C2 = table03(:, 2);

MapList02 = cell(length(C1), 2);

[~, Lia1] = ismember(C1, S2_ID);
[~, Lia2] = ismember(C2, S3_ID);

MapList02(:, 1) = info2.locus_tag(Lia1);
MapList02(:, 2) = info2.protein(Lia1);
MapList02(:, 3) = info2.gene(Lia1);

MapList02(:, 4) = info3.locus_tag(Lia2);
MapList02(:, 5) = info3.protein(Lia2);
MapList02(:, 6) = info3.gene(Lia2);
MapList02(:, 7)=table03(:,11);%E value
%% replace gene id
% Initialize new gene variable

% Extract and map relevant columns from MapList01
D01 = MapList01(:, 1);
D02 = MapList01(:, 4);
L1 = ismember(D01, Origin_gene);
A01 = [D01(L1) D02(L1)];

% Extract and map relevant columns from MapList02
E01 = MapList02(:, 1);
E02 = MapList02(:, 4);
L2 = ismember(E01, Origin_gene);
A11 =[ E01(L2) E02(L2)];
%% Find not map genes
% Find positions where L1 and L2 are zero
[~, L1] = ismember(Origin_gene,D01);
[~, L2] = ismember(Origin_gene,E01);

zeroPositions01 = find(L1 == 0);
zeroPositions02 = find(L2 == 0);

C = intersect(zeroPositions01,zeroPositions02);

% Extract relevant genes and proteins
RestGene01 = [model.proteins(C), model.genes(C)];

% Find matching genes in info3
[~, L3] = ismember(lower(info3.gene), lower(RestGene01(:, 1)));
matchingInfo3Idx = L3 ~= 0;
temp01 = [info3.gene(matchingInfo3Idx), info3.locus_tag(matchingInfo3Idx)];

% Find matching genes in RestGene01
[~, L4] = ismember(lower(temp01(:, 1)), lower(RestGene01(:, 1)));
temp02 = [RestGene01(L4, 2), info3.locus_tag(matchingInfo3Idx)];

% Find genes in RestGene01 that are not in temp02
RestGene02 = setdiff(RestGene01(:, 2), RestGene01(L4, 2));
RestGene02Map={'ECD_00709' 'K5T46_RS06905';'ECD_02645' 'K5T46_RS15505';'ECD_03483' 'K5T46_RS11605'};% 'ECD_RS18195' 'K5T46_RS11605'
% Combine the matrices into one matrix
ALL = [A01 ; A11 ; temp02;RestGene02Map];
ALL01=[ALL(:)];
m = size(ALL01, 1);
ALL = reshape(ALL01, [m/2, 2]);
[~, index01, indices01] =unique(ALL(:,1));
[~, index02, indices02] =unique(ALL(:,2));
A=[indices01 indices02];
[B,IA,IC]= unique(A, 'rows', 'stable');
newALL=ALL(IA,:);

%% Find unique elements in the first column and their indices
[uniqueElements, ~, indices] = unique(newALL(:, 1));

% Count occurrences of each unique element
elementCounts = histcounts(indices, 'BinMethod', 'integers')';

% Separate indices for single and multiple occurrences
singleOccurrenceIdx = find(elementCounts == 1);
multiOccurrenceIdx = find(elementCounts > 1);

% Retrieve rows with unique single occurrences
SingleIdx = ismember(indices, singleOccurrenceIdx);
GeneSingle = newALL(SingleIdx, :);

% Retrieve rows with multiple occurrences
MultiIdx = ismember(indices, multiOccurrenceIdx);
GeneMulti = newALL(MultiIdx, :);

%% replace gene id
[singlegeneslist, ~, ~] =unique(GeneSingle(:,1));
[multigeneslist, ~, ~] =unique(GeneMulti(:,1));

model=modelOut;
% Map new gene names from GeneSingle to NG
[LIA,LOCB]= ismember(Origin_gene,GeneSingle(:,1)) ;
NG=Origin_gene;
NG((LIA))=GeneSingle(LOCB(LIA),2);
model.genes=NG;
model.geneNames=NG;


model = creategrRulesField(model, model.rxns);
model= updateGenes(model);
model = generateRules(model);

%% update grRules
[results, ~] = findRxnsFromGenes(modelOut,singlegeneslist);
mergedResults=MergeStruct(results);
GeneRuleOld=mergedResults(:,5);     
GeneRuleNew=GeneRuleOld;

for i=1:length(GeneSingle)        
    GeneRuleNew=strrep(GeneRuleNew,GeneSingle{i,1},GeneSingle{i,2});    
end

for i=1:length(GeneSingle)  
     model = changeGeneAssociation(model, mergedResults(i,1),GeneRuleNew(i,1));     
end

for i=1:length(multigeneslist)
tempgenes=multigeneslist(i);
[results, ~] = findRxnsFromGenes(model,tempgenes);
mergedResults=MergeStruct(results);
Idx = ismember(GeneMulti(:,1), tempgenes);
ReplaceGenes=GeneMulti(Idx,2);
 for m=1:size(mergedResults, 1)% rxns numbers
     GeneRuleOld=mergedResults(m,5);
     rxnName=mergedResults(m,1);
     GeneRuleNew={};
     for n=1:length(ReplaceGenes)
         GeneRuleNew(n)=strrep(GeneRuleOld,tempgenes,ReplaceGenes(n));
     end
    combinedStr = strjoin(GeneRuleNew, ' or ');
    ReplaceStr = ['(' combinedStr ')'];
    model = changeGeneAssociation(model,rxnName,ReplaceStr);
 end
end

model = creategrRulesField(model, model.rxns);
[model] = updateGenes(model);
model = generateRules(model);
model.geneNames=model.genes;

%% define compartment
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

%% update other names

X=model.geneNames;
[~, Loc1] = ismember(X,info3.locus_tag);
newmodel = rmfield(model,{'geneisncbigiID','geneisrefseq_old_locus_tagID','proteinisrefseq_old_locus_tagID','proteinisncbigiID'});
newmodel.geneisrefseq_locus_tagID=model.genes;
newmodel.geneisrefseq_nameID=info3.gene(Loc1);
newmodel.proteinisrefseq_nameID=info3.protein_id(Loc1);
newmodel.proteinisrefseq_locus_tagID=model.geneNames;
newmodel.proteins=info3.protein(Loc1);
newmodel.modelID='NewBL21DE3';


outmodel = writeCbModel(newmodel, 'format','mat','fileName', 'NewBL21DE3.mat');



% sbmlModel01 = writeSBML(newmodel, fileName, compSymbolList, compNameList);

