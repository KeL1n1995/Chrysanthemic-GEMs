

function newmodel= ReviseMetsFormula(model,metsID, newformula)
[~,Met_Index] = ismember(metsID,model.mets);
model.metFormulas(Met_Index)=newformula;
newmodel=model;
end