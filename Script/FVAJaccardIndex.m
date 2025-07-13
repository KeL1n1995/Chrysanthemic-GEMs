
%%  Data source: Fig.2C Sorting results of FVA Jaccard index between BL21_WT and BL21_CS3
% Optional:
% [minFlux1, maxFlux1] = fluxVariability(model1, 90);
% [minFlux2, maxFlux2] = fluxVariability(model2, 90);


sample1=X1_1;
sample2=X1_2;

minFlux1=min(sample1,[],2);
maxFlux1=max(sample1,[],2);

minFlux2=min(sample2,[],2);
maxFlux2=max(sample2,[],2);

J = fvaJaccardIndex([minFlux1, minFlux2], [maxFlux1, maxFlux2]);
fprintf('Mean Jaccard index = %.4f.\n', mean(J));



%% 
% To visualise the FVA results, we plot the flux ranges as errorbars, with reactions 
% sorted by the Jaccard index.

highContrastOrange = [0.8500, 0.3250, 0.0980];
highContrastPurple = [0.4940, 0.1840, 0.5560];
highContrastBlue = [0, 0.4470, 0.7410];
highContrastGreen = [0.4660, 0.6740, 0.1880];
highContrastRed = [0.6350, 0.0780, 0.1840];
highContrastCyan = [0.3010, 0.7450, 0.9330];


clr1 = highContrastPurple;
clr2 =highContrastBlue;
clr3 =highContrastOrange;


E = [(maxFlux1 - minFlux1)/2 (maxFlux2 - minFlux2)/2];
Y = [minFlux1 minFlux2] + E;
X = [(1:length(Y)) - 0.1; (1:length(Y)) + 0.1]';
[~, xj] = sort(J);


f1 = figure;
    %errorbar(X, Y(xj, :), E(xj, :), 'linestyle', 'none', 'linewidth', 2);
hold on
    errorbar(X(:,1), Y(xj, 1), E(xj, 1), 'linestyle', 'none', 'linewidth', 1,'Color',clr1);
    errorbar(X(:,2), Y(xj, 2), E(xj, 2), 'linestyle', 'none', 'linewidth',1,'Color',clr2);
set(gca, 'xlim', [0, length(Y) + 1])
 
xlabel('Reaction')
ylabel('Flux range (mmol/gDW/h)')
ylim([-1000,1000])
yyaxis right

plot(J(xj),'linewidth', 1,'Color',clr3)
set(gca, 'YColor', clr3);

legend('WT', 'CS3', 'Jaccard','location', 'northoutside', ...
       'orientation', 'horizontal')
ylabel('Jaccard index')
hold off

hist(J);

Sample1=[min(sample1,[],2) mean(sample1,2) max(sample1,[],2)];
Sample2=[min(sample2,[],2) mean(sample2,2) max(sample2,[],2)];


hold on
pointSize = 4;
for i = 1:2000
        % Generate random x-coordinates around the center of the bar for scatter plot
        x = i ;
        scatter(x, Sample1(i,:),pointSize, 'filled', 'MarkerEdgeColor', clr_6, 'MarkerFaceColor', 'r');
end
for i = 1:2000
        % Generate random x-coordinates around the center of the bar for scatter plot
        x = i ;
        scatter(x, Sample2(i,:),pointSize, 'filled', 'MarkerEdgeColor', clr_55, 'MarkerFaceColor', 'r');
end
hold off


pointSize = 4;
scatter(Sample1(:,2),Sample2(:,2),pointSize, 'filled', 'MarkerEdgeColor', clr_6);




%% Data matrices
data1 = sample1;
data2 = sample2;

% Calculate statistics for data1
minValues1 = min(data1, [], 2);       % Minimum values
maxValues1 = max(data1, [], 2);       % Maximum values
medianValues1 = median(data1, 2);     % Median values
Q1Values1 = prctile(data1, 25, 2);    % First quartile (Q1)
Q3Values1 = prctile(data1, 75, 2);    % Third quartile (Q3)
meanValues1 = mean(data1, 2);         % Mean values

% Calculate statistics for data2
minValues2 = min(data2, [], 2);       % Minimum values
maxValues2 = max(data2, [], 2);       % Maximum values
medianValues2 = median(data2, 2);     % Median values
Q1Values2 = prctile(data2, 25, 2);    % First quartile (Q1)
Q3Values2 = prctile(data2, 75, 2);    % Third quartile (Q3)
meanValues2 = mean(data2, 2);         % Mean values

% --- Plotting ---
figure;
pointSize = 6; % Marker size for scatter plots

% Subplot 1: Min Values
subplot(3, 2, 1);
scatter(minValues1, minValues2, pointSize, 'filled', 'MarkerEdgeColor', highContrastOrange, 'MarkerFaceColor', highContrastOrange);
hold on;
plot([min(minValues1), max(minValues1)], [min(minValues1), max(minValues1)], 'k'); % Add y=x line
p = polyfit(minValues1, minValues2, 1); % Fit linear regression
yfit = polyval(p, minValues1);
plot(minValues1, yfit, 'r-'); % Plot regression line
title('Min Values');
xlabel('Min Values in WT_Sampling');
ylabel('Min Values 2');
hold off;

% Subplot 2: Max Values
subplot(3, 2, 2);
scatter(maxValues1, maxValues2, pointSize, 'filled', 'MarkerEdgeColor', highContrastPurple, 'MarkerFaceColor', highContrastPurple);
hold on;
plot([min(maxValues1), max(maxValues1)], [min(maxValues1), max(maxValues1)], 'k'); % Add y=x line
p = polyfit(maxValues1, maxValues2, 1); % Fit linear regression
yfit = polyval(p, maxValues1);
plot(maxValues1, yfit, 'r-'); % Plot regression line
title('Max Values');
xlabel('Max Values 1');
ylabel('Max Values 2');
hold off;

% Subplot 3: Median Values
subplot(3, 2, 3);
scatter(medianValues1, medianValues2, pointSize, 'filled', 'MarkerEdgeColor', highContrastBlue, 'MarkerFaceColor', highContrastBlue);
hold on;
plot([min(medianValues1), max(medianValues1)], [min(medianValues1), max(medianValues1)], 'k'); % Add y=x line
p = polyfit(medianValues1, medianValues2, 1); % Fit linear regression
yfit = polyval(p, medianValues1);
plot(medianValues1, yfit, 'r-'); % Plot regression line
title('Median Values');
xlabel('Median Values 1');
ylabel('Median Values 2');
hold off;

% Subplot 4: Q1 Values
subplot(3, 2, 4);
scatter(Q1Values1, Q1Values2, pointSize, 'filled', 'MarkerEdgeColor', highContrastGreen, 'MarkerFaceColor', highContrastGreen);
hold on;
plot([min(Q1Values1), max(Q1Values1)], [min(Q1Values1), max(Q1Values1)], 'k'); % Add y=x line
p = polyfit(Q1Values1, Q1Values2, 1); % Fit linear regression
yfit = polyval(p, Q1Values1);
plot(Q1Values1, yfit, 'r-'); % Plot regression line
title('Q1 Values');
xlabel('Q1 Values 1');
ylabel('Q1 Values 2');
hold off;

% Subplot 5: Q3 Values
subplot(3, 2, 5);
scatter(Q3Values1, Q3Values2, pointSize, 'filled', 'MarkerEdgeColor', highContrastRed, 'MarkerFaceColor', highContrastRed);
hold on;
plot([min(Q3Values1), max(Q3Values1)], [min(Q3Values1), max(Q3Values1)], 'k'); % Add y=x line
p = polyfit(Q3Values1, Q3Values2, 1); % Fit linear regression
yfit = polyval(p, Q3Values1);
plot(Q3Values1, yfit, 'r-'); % Plot regression line
title('Q3 Values');
xlabel('Q3 Values 1');
ylabel('Q3 Values 2');
hold off;

% Combine statistics into new matrices
BP_data1 = [minValues1 maxValues1 medianValues1 Q1Values1 Q3Values1];
BP_data2 = [minValues2 maxValues2 medianValues2 Q1Values2 Q3Values2];