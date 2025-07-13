
%% Data source: Fig.2D Volcano plot of reaction fluxes from flux sampling in the WT and CS3 models. 

rxn_name=Model.rxns;

data1=X1_1;% Variable X1_1 can be obtained by running the script SamplingAndFVA
data2=X1_2;% Variable X1_2 can be obtained by running the script SamplingAndFVA

meanValues1 = mean(data1,  2);
meanValues2 = mean(data2,  2);

indices1 = find(abs(meanValues1) == 0);
meanValues1(indices1) = 10^-8;
indices2 = find(abs(meanValues2) ==0);
meanValues2(indices2) = 10^-8;

FC=meanValues2./meanValues1;% cs3/wt

[rowNum, colNum] = size(data1);
pValues = zeros(rowNum, 1); % 
% random sampling
sampleSize = 300;
indices = randperm(size(data1, 1), sampleSize);
X_sample = data1(:, indices);
Y_sample = data2(:, indices);

for i = 1:rowNum
[~, pValues(i,1)]= ttest2(X_sample(i, :), Y_sample(i, :));
end

pValues(pValues==0)=10^-8;

indices3 = find(~isnan(FC)&(FC>0));

indices4 = find(~isnan(pValues)&(pValues>0));

idx=intersect(indices3,indices4);

LG_pValues=-log10(pValues(idx));
LG_FC=log2(FC(idx));

rxn_name=Model.rxns(idx);

%%
genes = rxn_name;
log2FC = LG_FC;% Log2 fold change
%  -log10 p-value
negLog10P = LG_pValues;

colors = zeros(length(log2FC), 3); 
rgb1 = [0, 0.7529, 0];
rgb2 = [0.949, 0.714, 0.506];

colors = repmat([0.5, 0.5, 0.5], length(log2FC), 1);
colors(log2FC > 0, :) = repmat(rgb2, sum(log2FC > 0), 1); % 
colors(log2FC <= 0, :) = repmat(rgb1, sum(log2FC <= 0), 1); % 


figure;
scatter(log2FC, negLog10P, 25, 'MarkerEdgeColor', 'flat', 'MarkerFaceColor', 'none', 'CData', colors); % 
hold on;


Sel_idx=[710;711;960;2285;2311];


for n=1:length(Sel_idx)
i=Sel_idx(n)

scatter(log2FC(i), negLog10P(i), 25, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r', 'LineWidth', 1.5);
% text(log2FC(i), negLog10P(i), genes{i}, 'HorizontalAlignment', 'right', ...
%     'VerticalAlignment', 'bottom', 'FontSize', 10, 'Color', 'black');
end

xlabel('Log_2 Fold Change');

xlim([-8, 8]);

ylabel('-Log_{10} P-value');
% title('Volcano Plot');
grid on;
set(gca, 'FontSize', 12);


pThreshold = 10^-8;
fcThreshold = 1;
line(xlim, [-log10(pThreshold), -log10(pThreshold)], 'Color', 'k', 'LineStyle', '--', 'DisplayName', 'p = 10^-8');
line([fcThreshold, fcThreshold], ylim, 'Color', rgb2, 'LineStyle', '--', 'DisplayName', 'Log2FC = 1');
line([-fcThreshold, -fcThreshold], ylim, 'Color', rgb1, 'LineStyle', '--', 'DisplayName', 'Log2FC = -1');

set(groot, 'DefaultAxesFontName', 'Arial'); % 
set(groot, 'DefaultTextFontName', 'Arial'); % 

