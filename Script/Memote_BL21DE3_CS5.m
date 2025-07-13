
%% Memote test
model=BL21_CS5;


model.metFormulas = strrep(model.metFormulas,'P0','P');% Replace 'P0' in phosphorus-containing metabolite formulas
rxnsToIgnore = ones(length(model.rxns),1); % Do not consider exchange, source, sink, and biomass reactions since they are inherently imbalanced
rxnsToIgnore((findExcRxns(model))) = 0;
rxnsToIgnore((model.c==1)) = 0;
model.SIntRxnBool = logical(rxnsToIgnore);
[massImbalance,imBalancedMass,imBalancedCharge,imBalancedRxnBool,Elements,missingFormulaeBool,balancedMetBool] = checkMassChargeBalance(model,-1);
imbalancedRxnsMass = setdiff(find(~cellfun(@isempty,imBalancedMass)),find(rxnsToIgnore == 0));


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


outmodel = writeCbModel(model, 'format','mat', 'fileName', 'NewBL21DE3_CS5.mat');


 %% Optional : get memote report
% sbmlModel = writeSBML(model, 'NewBL21DE3_CS5');
% command = 'memote report snapshot --filename "Memote_BL21DE3_CS5.html" NewBL21DE3_CS5.xml'; % 
% status = system(command);
