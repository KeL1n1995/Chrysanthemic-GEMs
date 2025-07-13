
%% Predicted enzyme reaction flux as a function of substrate concentration and gene expression (Fig.4D)

%% ADH2
% --- 1. Define Enzyme Kinetic Parameters ---
kcat_ADH2 = 0.75;      % Kcat
Km_ADH2 = 236.0e-03;      % Km 

% --- 2. Define Approximate Conversion Factor from Gene Expression to Enzyme Concentration ---

% This is a hypothetical value used to convert Enzyme Gene Expression(EXP)
% into molar enzyme concentration (mM).
% In real applications, this requires experimental data (e.g., proteomics) for calibration.

exp_mRNA1=40836.24;% from experimental data
alpha_ADH2=f1/kcat_ADH2/exp_mRNA1;% Hypothetical conversion factor (M / EXP unit) ; 
%Parameter "f1" is obtained by running script "SimulateEnhanceGeneExpression"
alpha=alpha_ADH2;

% --- 3. Define the Range for Variables ---
% Substrate Concentration [S]
S_min = 0;
S_max =2*Km_ADH2;   % Maximum substrate concentration (0.5 mM)
num_S_points = 50; % Number of points on the S-axis
S = linspace(S_min, S_max, num_S_points);

% Enzyme Gene Expression
EXP_min = 0;
EXP_max = 100000; % Maximum TPM value
num_EXP_points = 50; % Number of points on the TPM-axis
EXP = linspace(EXP_min, EXP_max, num_EXP_points);

% --- 4. Create Mesh Grid Data ---
% This creates coordinate pairs for all combinations of S and TPM
[S_mesh1, EXP_mesh1] = meshgrid(S, EXP);

% --- 5. Calculate Corresponding Enzyme Concentration [E]_T ---
% Convert TPM to enzyme concentration using the assumed alpha
Et_mesh1 = alpha * EXP_mesh1;

% --- 6. Calculate Reaction Flux (V) ---

V_flux1 = (kcat_ADH2 .* Et_mesh1 .* S_mesh1) ./ (Km_ADH2 + S_mesh1);

%% Plot the 3D Figure ---

figure;

h = surf(S_mesh1, EXP_mesh1, V_flux1);

title('ADH2 Enzyme');
xlabel('Substrate Concentration [S] (mM)');
ylabel('Gene Expression');
zlabel('Reaction Flux V (mM\cdot{}s^{-1})', 'Interpreter', 'tex');
zlabel_handle = get(gca,'ZLabel'); % Get its handle
set(zlabel_handle, 'FontName', 'Arial', 'FontSize', 12); % Set its font properties
% Adjust view angle and color settings
view(-37.5, 30)
colorbar; % Display the color bar, indicating Z-axis values
colormap('cool'); % Set the colormap
shading interp; % Smooth colors

% Make the surface translucent
alpha_value = 0.6; 
set(h, 'FaceAlpha', alpha_value);

% Add grid lines for better visualization
grid on;

% --- Set Font Properties ---
set(gca, 'FontName', 'Arial', 'FontSize', 12); % For axis labels and ticks
set(get(gca,'Title'), 'FontName', 'Arial', 'FontSize', 12); % For title
set(get(gca,'XLabel'), 'FontName', 'Arial', 'FontSize', 12); % For x-label
set(get(gca,'YLabel'), 'FontName', 'Arial', 'FontSize', 12); % For y-label
set(get(gca,'ZLabel'), 'FontName', 'Arial', 'FontSize', 12); % For z-label

% For the legend (if it exists)
legend_handle = findobj(gcf, 'Type', 'Legend');
if ~isempty(legend_handle)
    set(legend_handle, 'FontName', 'Arial', 'FontSize', 12);
end

% % --- Save the figure (as high-res PNG and TIFF) ---
resolution_dpi = 1200; % Use a high DPI for clarity

% Save as PNG
png_filename = 'ADH2_Flux_3D_Plot.png';
print(gcf, png_filename, '-dpng', ['-r', num2str(resolution_dpi)]);
disp(['Plot saved as: ' png_filename]);

% % Save as TIFF
% tiff_filename = 'ADH2_Flux_3D_Plot.tif';
% print(gcf, tiff_filename, '-dtiff', ['-r', num2str(resolution_dpi)]);
% disp(['Plot saved as: ' tiff_filename]);

% filename = 'ADH2_Flux_3D_Plot.pdf'; % Choose a descriptive name and .pdf extension
% % Save the current figure (gcf)
% % '-vector' ensures it's saved as a vector graphic if the format supports it
% % '-painters' is often recommended for vector graphics, especially with transparency
% print(gcf, filename, '-dpdf', '-vector', '-vector');
% disp(['Plot saved as: ' filename]);

%% ALDH1
% --- 1. Define Enzyme Kinetic Parameters ---
kcat_ALDH1 = 0.11;      % Kcat 
Km_ALDH1 = 4.4e-03;      % Km 

% --- 2. Define Approximate Conversion Factor from Gene Expression to Enzyme Concentration ---

% This is a hypothetical value used to convert Enzyme Gene Expression(EXP)
% into molar enzyme concentration (mM).
% In real applications, this requires experimental data (e.g., proteomics) for calibration.

exp_mRNA2=48677.04218;% from experimental data
alpha_ALDH1=f2/kcat_ALDH1/exp_mRNA2;% Hypothetical conversion factor (M / EXP unit)
alpha=alpha_ALDH1;

% --- 3. Define the Range for Variables ---
% Substrate Concentration [S]
S_min = 0;
S_max =2*Km_ALDH1;   % Maximum substrate concentration (0.5 mM)
num_S_points = 50; % Number of points on the S-axis
S = linspace(S_min, S_max, num_S_points);

% Enzyme Gene Expression
EXP_min = 0;
EXP_max = 100000; % Maximum TPM value
num_EXP_points = 50; % Number of points on the TPM-axis
EXP = linspace(EXP_min, EXP_max, num_EXP_points);

% --- 4. Create Mesh Grid Data ---
% This creates coordinate pairs for all combinations of S and EXP
[S_mesh2, EXP_mesh2] = meshgrid(S, EXP);

% --- 5. Calculate Corresponding Enzyme Concentration [E]_T ---
% Convert TPM to enzyme concentration using the assumed alpha
Et_mesh2 = alpha * EXP_mesh2;

% --- 6. Calculate Reaction Flux (V) ---

V_flux2 = (kcat_ALDH1 .* Et_mesh2 .* S_mesh2) ./ (Km_ALDH1 + S_mesh2);

%%  Plot the 3D Figure ---
resolution_dpi = 1200;
figure;

h = surf(S_mesh2, EXP_mesh2, V_flux2);

title('ALDH1 Enzyme');
xlabel('Substrate Concentration [S] (mM)');
ylabel('Gene Expression');
zlabel('Reaction Flux V (mM\cdot{}s^{-1})', 'Interpreter', 'tex');
zlabel_handle = get(gca,'ZLabel'); % Get its handle
set(zlabel_handle, 'FontName', 'Arial', 'FontSize', 12); % Set its font properties
% Adjust view angle and color settings
% view(3); % 3D view

view(-37.5, 30)
colorbar; % Display the color bar, indicating Z-axis values
colormap('cool'); % Set the colormap
shading interp; % Smooth colors

% Make the surface translucent
alpha_value = 0.6; % You can adjust this value (0 for fully transparent, 1 for fully opaque)
set(h, 'FaceAlpha', alpha_value); % Apply transparency to the surface object

% Add grid lines for better visualization
grid on;

% --- Set Font Properties ---
set(gca, 'FontName', 'Arial', 'FontSize', 12); % For axis labels and ticks
set(get(gca,'Title'), 'FontName', 'Arial', 'FontSize', 12); % For title
set(get(gca,'XLabel'), 'FontName', 'Arial', 'FontSize', 12); % For x-label
set(get(gca,'YLabel'), 'FontName', 'Arial', 'FontSize', 12); % For y-label
set(get(gca,'ZLabel'), 'FontName', 'Arial', 'FontSize', 12); % For z-label

% For the legend (if it exists)
legend_handle = findobj(gcf, 'Type', 'Legend');
if ~isempty(legend_handle)
    set(legend_handle, 'FontName', 'Arial', 'FontSize', 12);
end

% filename = 'ALDH1_Flux_3D_Plot.pdf'; % Choose a descriptive name and .pdf extension
% 
% % Save the current figure (gcf)
% % '-vector' ensures it's saved as a vector graphic if the format supports it
% % '-painters' is often recommended for vector graphics, especially with transparency
% print(gcf, filename, '-dpdf', '-vector', '-vector');
% disp(['Plot saved as: ' filename]);

% % --- Save as TIFF ---
% tiff_filename = 'ALDH1_Flux_3D_Plot.tif';
% % '-dtiff' for TIFF format
% % '-r<DPI>' to specify resolution
% print(gcf, tiff_filename, '-dtiff', ['-r', num2str(resolution_dpi)]); 
% disp(['Plot saved as: ' tiff_filename]);

% --- Save as PNG ---
png_filename = 'ALDH1_Flux_3D_Plot.png';
% '-dpng' for PNG format
print(gcf, png_filename, '-dpng', ['-r', num2str(resolution_dpi)]);
disp(['Plot saved as: ' png_filename]);



