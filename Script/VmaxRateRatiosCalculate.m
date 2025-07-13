%%  Data source: Fig. S15. The catalytic rate ratios of different enzymes interacting with the DMAPP.

Table_UniKP=readtable('Database/UniKP_Result.xlsx');
load('TranscriptomicData.mat', 'EXPO')
Kcat=Table_UniKP.Kcat;
Km=Table_UniKP.Km;
Name=Table_UniKP.ID;

S = logspace(-3, 3, 7);


TPM=EXPO.CS6;
GeneID=EXPO.Gene;
[C,IA,IC]=intersect(GeneID,Name);

E=mean(TPM(IA,:),2);

Km_1=Km(IC);
Kcat_1=Kcat(IC);
Ratio=[];


for x=1:length(Kcat_1)    
    for y=1:7
Test=(E(x)*Kcat_1(x))/(Km_1(x)+S(y));    
Total = sum((E.*Kcat_1 )./ (Km_1 + S(y)));
Ratio(x,y)=Test/Total;
    end
end

Name(IC)

