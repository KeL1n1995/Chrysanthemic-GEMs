% Example data
categories = {'WT-Geraniol', 'CS3-Geraniol', 'WT-Farnesol', 'CS3-Farnesol','WT-Farnesal', 'CS3-Farnesal'};
data=[X1_1(idx1, :);X1_2(idx1, :);X1_1(idx2, :);X1_2(idx2, :);X1_1(idx3, :);X1_2(idx3, :)];
summaryData=mean(data,2);
individualData=data;

% Number of categories
numCategories = length(categories);

% Create a new figure
figure;

% Create the bar plot
% b = bar(summaryData, 'FaceColor', 'flat');

% Hold the plot to overlay the scatter plot
hold on;
   
pointSize = 1; 

clr_6 = [5,113,176]/255; 
clr_55 = [202,0,32]/255; 

    % Plot each column as a separate scatter plot
    for i = 1:2:6
        % Generate random x-coordinates around the center of the bar for scatter plot
        x = i + 0.1 * randn(size(individualData(1,:)));
        scatter(x, individualData(i,:),pointSize, 'filled', 'MarkerEdgeColor', clr_6, 'MarkerFaceColor', 'r');
    end

   for i = 2:2:6
        % Generate random x-coordinates around the center of the bar for scatter plot
        x = i + 0.1 * randn(size(individualData(1,:)));
        scatter(x, individualData(i,:),pointSize, 'filled', 'MarkerEdgeColor',clr_55, 'MarkerFaceColor', 'r');
    end

    % Customize plot
    
    xx=1:6;
    set(gca, 'XTick', xx, 'XTickLabel', categories);

    xlabel('Compound');
    ylabel('Flux (mmol/gCDW/h)');
    
%     title('Combined Bar and Scatter Plot');
%     legend('Summary Data', 'Individual Data');

    % Optional: Add grid lines
    
    grid on;
    
    % Release the hold on the plot
    
    hold off;
    
    hold on;
     for i = 1:2746
    scatter(X1_1(i, :),X1_2(i, :),pointSize, 'filled', 'MarkerEdgeColor',clr_55, 'MarkerFaceColor', 'r');
     end
     hold off;
    
     