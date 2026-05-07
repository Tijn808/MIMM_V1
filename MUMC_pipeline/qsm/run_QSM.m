% QSM reconstruction for MUMC ME-GRE data using SEPIA
%
% Acquisition parameters (fixed for MUMC protocol):
%   Sequence  : 3D ME-GRE (monopolar readout)
%   Scanner   : 3T Siemens
%   TE        : 6, 12, 18, 24, 30 ms
%   TR        : 35 ms
%   Voxel     : 1 x 1 x 1 mm isotropic
%   Matrix    : 224 x 224 x 154
%   B0_dir    : [0, 0, 1]  (standard axial acquisition)
%
% Requirements:
%   - SEPIA toolbox (https://github.com/kschan0214/sepia)
%   - Input: 4D magnitude NIfTI  [224 x 224 x 154 x 5]
%            4D phase NIfTI      [224 x 224 x 154 x 5]
%
% Output:
%   - QSM map     : NIfTI [224 x 224 x 154], units ppm
%   - Brain mask  : NIfTI [224 x 224 x 154], binary

%% --- Acquisition parameters (do not change for MUMC data) ---

TE       = [6, 12, 18, 24, 30] * 1e-3;   % echo times in seconds
B0       = 3;                              % field strength in Tesla
B0_dir   = [0, 0, 1];                     % B0 direction (axial)
voxel    = [1, 1, 1];                     % voxel size in mm

%% --- Input/output paths (adapt per subject) ---

% Set subject directory before running
subj_dir    = '';           % e.g. '/data/MUMC/sub-01/'
magnitude   = fullfile(subj_dir, 'magnitude.nii.gz');
phase       = fullfile(subj_dir, 'phase.nii.gz');
output_dir  = fullfile(subj_dir, 'qsm/');

if ~exist(output_dir, 'dir'); mkdir(output_dir); end

%% --- Load data ---

mag_info  = niftiinfo(magnitude);
mag_data  = niftiread(mag_info);    % [224 x 224 x 154 x 5]

pha_info  = niftiinfo(phase);
pha_data  = niftiread(pha_info);    % [224 x 224 x 154 x 5]
% Scale phase to radians if stored as integer (check header units)
% pha_data = double(pha_data) * pha_info.MultiplicativeScaling;

%% --- Step 1: Brain mask ---
% Use the first echo magnitude for brain extraction.
% FSL BET is called externally; the mask is loaded here.
% Run in terminal: bet magnitude_echo1 brain -m -f 0.4

brain_mask_file = fullfile(output_dir, 'brain_mask.nii.gz');
if ~exist(brain_mask_file, 'file')
    error('Brain mask not found. Run FSL BET first: bet %s %s -m -f 0.4', ...
          magnitude, fullfile(output_dir, 'brain'));
end
mask_info = niftiinfo(brain_mask_file);
Brain_Mask = logical(niftiread(mask_info));

%% --- Step 2: Total field estimation ---
% Multi-echo field map estimation via SEPIA
% Fits phase evolution across TEs to estimate field map and R2*

sepia_header.TE      = TE;
sepia_header.B0      = B0;
sepia_header.B0_dir  = B0_dir(:);
sepia_header.voxelSize = voxel;

% SEPIA field map estimation (echo combination + unwrapping)
% Second output is R2* map (s^-1) — saved for chi-separation
[field_map, R2star, ~, ~] = estimateTotalField(pha_data, mag_data, Brain_Mask, ...
    sepia_header, ...
    'method', 'MEDI nonlinear fit', ...
    'Unwrap',  'Laplacian');

%% --- Step 3: Background field removal ---
% Remove contributions from outside the brain (air, skull, shimming)
% VSHARP is robust for brain QSM

[local_field, ~] = removeBFwithVSHARP(field_map, Brain_Mask, voxel, ...
    'radius', 1:10);    % spherical kernel radii 1–10 mm

%% --- Step 4: Dipole inversion ---
% Reconstruct QSM from local field map via MEDI

QSM = MEDI_L1('field',   local_field, ...
               'mask',    Brain_Mask, ...
               'voxel',   voxel, ...
               'B0_dir',  B0_dir, ...
               'B0',      B0, ...
               'magnitude', mag_data(:,:,:,1));

%% --- Save outputs ---

% QSM map (ppm)
qsm_info            = mag_info;
qsm_info.Datatype   = 'double';
qsm_info.ImageSize  = size(QSM);
qsm_info.PixelDimensions = voxel;
niftiwrite(QSM, fullfile(output_dir, 'QSM.nii.gz'), qsm_info, 'Compressed', true);

% R2* map (s^-1) — needed for chi-separation
r2s_info          = qsm_info;
niftiwrite(R2star, fullfile(output_dir, 'R2star.nii.gz'), r2s_info, 'Compressed', true);

% Brain mask
mask_out            = qsm_info;
mask_out.Datatype   = 'uint8';
niftiwrite(uint8(Brain_Mask), fullfile(output_dir, 'brain_mask.nii.gz'), mask_out, 'Compressed', true);

disp('QSM reconstruction complete.')
disp(['Output written to: ' output_dir])
