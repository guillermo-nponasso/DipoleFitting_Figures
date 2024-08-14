%% clear environment

restoredefaultpath;
clear all;
addpath(fullfile('../Engine'));

%% global options

subject='117122';

%% load planes

p1 = load(fullfile('data',sprintf("%s_gm_plane.mat",subject)));   % grey matter cutoff plane
p2 = load(fullfile('data',sprintf("%s_bone_plane.mat",subject))); % skull cutoff plane

%% load stats

errormap = load(fullfile('data', sprintf("%s_interpolated_error.mat", subject)));
tSKIN = load(fullfile('data', sprintf("%s_thickness_SKIN.mat", subject)));
tBONE = load(fullfile('data', sprintf("%s_thickness_BONE.mat", subject)));
tCSF  = load(fullfile('data', sprintf("%s_thickness_CSF.mat", subject)));
tGM   = load(fullfile('data', sprintf("%s_thickness_GM.mat", subject)));

%% load GM shell

fprintf("Loading GM shell of subject: %s -- Please wait ..\n", subject);
GM = stlread(fullfile('data',sprintf("%s_gm_headreco.stl",subject)));
disp("GM shell loaded!");

GM_P = GM.Points;
GM_t = GM.ConnectivityList;

tricent = meshtricenter(GM_P,GM_t);
trinorm = meshnormals(GM_P,GM_t);

% check sizes of data
if length(tricent) == length(errormap.ErrGM)
    disp("Errormap data length check -- OK");
else
    disp("Errormap data length check -- FAIL");
end

%% load cluster data

CL = load(fullfile('data','clusters',sprintf("%s_cluster_4000.mat", subject)));

%% obtain offset and normal vectors of each plane

n1 = [0; -sin(p1.x_rad); cos(p1.x_rad)];
v1 = p1.cutoff * n1;

n2 = [0; -sin(p2.x_rad); cos(p2.x_rad)];
v2 = p2.cutoff * n2;

%% (optional) plot the restricted GM
% 
% below_tri_ix = (tricent-repmat(v1',length(tricent),1))*n1 < 0;
% vals = (tricent-repmat(v1',length(tricent),1))*n1;
% vals(below_tri_ix) = NaN;
% 
% h=figure;
% patch('faces', GM_t, 'vertices', GM_P, 'FaceVertexCData', vals, 'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceAlpha', 1.0);                   
% colormap('parula');
% colorbar;
% clim([min(vals) max(vals)]);
% axis 'equal';  axis 'tight';      
% xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');

%% filter by cutoff planes and compute averages

% error means ---
GM_tri_ix = (tricent-repmat(v1',length(tricent),1))*n1 >= 0;
r_error   = errormap.ErrGM(GM_tri_ix); % restricted error map

fprintf("===================================================\n");
fprintf("The average (interpolated) fitting error is: %.2d (mm)\n", mean(r_error));

CL_ix = (CL.cl_centers - repmat(v1', length(CL.cl_centers),1))*n1 >= 0;
r_error   = errormap.errors(CL_ix);

fprintf("The average fitting error is: %.2d (mm)\n", mean(r_error));
fprintf("===================================================\n");
% ---

% skin thickness ---
SKIN_ix  = (tSKIN.VNodes - repmat(v2', length(tSKIN.VNodes),1))*n2 >= 0;
tskin = tSKIN.VDist(SKIN_ix);

fprintf("The average skin thickness is: %.2d (mm)\n", mean(tskin));
fprintf("===================================================\n");
% ---

% bone thickness
BONE_ix  = (tBONE.VNodes - repmat(v2', length(tBONE.VNodes),1))*n2 >= 0;
tbone = tBONE.VDist(BONE_ix);

fprintf("The average bone thickness is: %.2d (mm)\n", mean(tbone));
fprintf("===================================================\n");
% ---

% CSF thickness
CSF_ix  = (tCSF.VNodes - repmat(v1', length(tCSF.VNodes),1))*n1 >= 0;
tcsf = tCSF.VDist(CSF_ix);
fprintf("The average CSF thickness is: %.2d (mm)\n", mean(tcsf));
fprintf("===================================================\n");
% ---

% GM thickness
GM_ix  = (tGM.VNodes - repmat(v1', length(tGM.VNodes),1))*n1 >= 0;
tgm = tGM.VDist(GM_ix);
fprintf("The average GM thickness is: %.2d (mm)\n", mean(tgm));
fprintf("===================================================\n");
% ---

%% create table and save
r_error(isnan(r_error))=[];
% Create table
subject_table = table(str2num(subject), mean(r_error), mean(tskin), mean(tbone), ...
    mean(tcsf), mean(tgm), ...
    'VariableNames', {'Subject', 'Mean_Fit_Error_mm', 'Mean_Skin_Thickness_mm', ...
    'Mean_Bone_Thickness_mm', 'Mean_CSF_Thickness_mm', 'Mean_GM_Thickness_mm'});

writetable(subject_table, fullfile('data',sprintf("%s_averages.csv", subject)));

disp("Written subject data.")

%% calculate the inter-tissue correlations

SKIN_ix_2 = SKIN_ix;
SKIN_ix_1 = (tSKIN.VNodes - repmat(v1', length(tSKIN.VNodes),1))*n1 >= 0;
SKIN_ix   = SKIN_ix_1 & SKIN_ix_2;

inflated_error = load(fullfile('data', sprintf('%s_interpolated_error_SKIN.mat', subject))).ErrSK;
inflated_GM    = load(fullfile('data', sprintf('%s_inflated_thickness_GM.mat', subject))).thickSK;
inflated_CSF   = load(fullfile('data', sprintf('%s_inflated_thickness_CSF.mat', subject))).thickSK;
inflated_BONE  = load(fullfile('data', sprintf('%s_inflated_thickness_BONE.mat', subject))).thickSK;

ierror = inflated_error(SKIN_ix);
igm    = inflated_GM(SKIN_ix);
icsf   = inflated_CSF(SKIN_ix);
ibone  = inflated_BONE(SKIN_ix);
iskin  = tSKIN.VDist(SKIN_ix);

data = [ierror igm icsf ibone iskin];
data(any(isnan(data), 2), :) = []; % remove nan from data

std_vec = std(data,0,1);
normalization_matrix = (std_vec' * std_vec).^(-1);
covariance_matrix = cov(data);
fprintf("Error - tissue thickness correlation matrix:\n");
correlation_matrix = normalization_matrix.*covariance_matrix

labels = ['fit_error', 'gm_thickness', 'csf_thickness', 'bone_thickness', 'skin_thickness'];
save(fullfile('data',sprintf('%s_error_tissue_cor_matrix.mat', subject)), 'correlation_matrix', 'labels');
