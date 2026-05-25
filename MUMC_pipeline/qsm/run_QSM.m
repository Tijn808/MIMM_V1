% QSM reconstruction for MUMC ME-GRE data using SEPIA built-in functions
%
% Acquisition parameters (fixed for MUMC protocol):
%   Sequence  : 3D ME-GRE (monopolar readout), Philips MR 7700
%   Scanner   : 3T
%   TE        : 6.001, 12, 18, 24, 30 ms
%   TR        : 35 ms
%   Voxel     : 1 x 1 x 1 mm isotropic
%   Matrix    : 224 x 224 x 154
%   B0_dir    : [0, 0, 1]
%
% Prerequisites:
%   Run MUMC_pipeline/preprocess/prepare_mgre.m first to produce:
%     magnitude.nii.gz  [224 x 224 x 154 x 5]
%     phase.nii.gz      [224 x 224 x 154 x 5]  (radians)
%
% Outputs (written to output_dir):
%   QSM.nii.gz          — susceptibility map (ppm)
%   R2star.nii.gz       — R2* map (s^-1)
%   brain_mask.nii.gz   — binary brain mask

%% --- Paths ---

sepia_dir  = '/home/tijn-saes/Documents/Internship/MIMM/sepia';
subj_dir   = '/home/tijn-saes/Documents/Internship/ME_GRE/';
output_dir = fullfile(subj_dir, 'qsm');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

%% --- Add SEPIA to path ---

addpath(genpath(sepia_dir));

%% --- Acquisition parameters ---

TE      = [0.006001, 0.012, 0.018, 0.024, 0.030];   % seconds
B0      = 3;                                          % Tesla
B0_dir  = [0, 0, 1];
voxel   = [1, 1, 1];                                  % mm
gamma   = 2 * pi * 42.577478518e6;                    % rad/s/T

%% --- Load data ---

disp('Loading magnitude and phase...')
mag_info  = niftiinfo(fullfile(subj_dir, 'magnitude.nii.gz'));
mag_data  = double(niftiread(mag_info));    % [224 x 224 x 154 x 5]
pha_data  = double(niftiread(fullfile(subj_dir, 'phase.nii.gz')));  % radians

%% --- Brain mask ---
% Use FSL BET on first echo magnitude. Saves brain_mask.nii.gz.

mask_file = fullfile(output_dir, 'brain_mask.nii.gz');
if ~exist(mask_file, 'file')
    disp('Running FSL BET for brain mask...')
    bet_input  = fullfile(output_dir, 'mag_e1');
    bet_output = fullfile(output_dir, 'brain');

    % Save echo 1 magnitude as temporary NIfTI for BET
    ref1             = mag_info;
    ref1.ImageSize   = size(mag_data(:,:,:,1));
    ref1.Datatype    = 'double';
    ref1.BitsPerPixel = 64;
    ref1.PixelDimensions = voxel;
    niftiwrite(mag_data(:,:,:,1), [bet_input '.nii.gz'], ref1, 'Compressed', true);

    system(sprintf('/home/tijn-saes/fsl/bin/bet %s %s -m -f 0.4', ...
        [bet_input '.nii.gz'], bet_output));
    % BET with -m writes <bet_output>_mask.nii.gz = brain_mask.nii.gz directly
end

Brain_Mask = logical(niftiread(mask_file));
fprintf('Brain mask: %d voxels (%.1f%%)\n', sum(Brain_Mask(:)), 100*mean(Brain_Mask(:)));

%% --- Multi-echo field map estimation ---
% Phase difference between adjacent echoes avoids large wraps.
% Weighted least-squares fit of phase vs TE gives field in rad/s.

disp('Estimating field map from multi-echo phase...')
matrixSize = size(Brain_Mask);

% Phase-difference field map: use angle(exp(i*dPhi)) per adjacent pair.
% Each pair spans only dTE=6 ms, so wrapping occurs only above ~83 Hz
% (>> typical brain susceptibility values), avoiding multi-wrap errors
% that corrupt an absolute-phase linear regression across all echoes.
dTE = TE(2) - TE(1);   % 6 ms, constant spacing

field_sum   = zeros(matrixSize);
weight_sum  = zeros(matrixSize);

for e = 1:length(TE)-1
    dphi    = angle(exp(1i * (pha_data(:,:,:,e+1) - pha_data(:,:,:,e))));
    field_e = dphi / dTE;                                    % rad/s
    w       = mag_data(:,:,:,e) .* mag_data(:,:,:,e+1);     % geometric mean weight
    field_sum  = field_sum  + w .* field_e;
    weight_sum = weight_sum + w;
end

field_map    = field_sum ./ max(weight_sum, eps);   % rad/s
field_map_hz = field_map / (2 * pi);                % Hz

%% --- Background field removal (VSHARP) ---
% Phase-difference field map is already unwrapped (dTE=6ms means wrapping
% only above ~83 Hz, well above brain tissue values). Skip unwrapping.

disp('Removing background field with VSHARP...')
[local_field, Brain_Mask_vsharp] = BKGRemovalVSHARP(field_map_hz, Brain_Mask, ...
    matrixSize, 'radius', 1:10);

%% --- Dipole inversion (iLSQR) ---

disp('Running dipole inversion (iLSQR)...')
QSM = qsmIterativeLSQR(local_field, Brain_Mask_vsharp, matrixSize, voxel, ...
    'b0dir', B0_dir);

% Convert from Hz to ppm
QSM = QSM / (B0 * gamma / (2 * pi)) * 1e6;

%% --- R2* map from magnitude decay ---

disp('Fitting R2*...')
% Vectorised R2* fit: log(S) = -R2* * TE + log(S0)
log_mag     = log(max(mag_data, 1e-6));        % avoid log(0)
log_mag_vec = reshape(log_mag, [], length(TE));
A_r2s       = [-TE(:), ones(length(TE), 1)];
coeff_r2s   = A_r2s \ log_mag_vec';            % [2 x N_voxels]
R2star      = reshape(coeff_r2s(1,:), matrixSize);
R2star      = max(R2star, 0) .* Brain_Mask;

%% --- Save outputs ---

% Bootstrap a clean 3D NIfTI header from a dummy write
dummy_file = fullfile(tempdir, 'mimm_ref3d.nii');
niftiwrite(single(QSM), dummy_file);
ref = niftiinfo(dummy_file);
delete(dummy_file);
ref.PixelDimensions = voxel;
ref.SpaceUnits      = 'Millimeter';

ref_dbl          = ref; ref_dbl.Datatype = 'double'; ref_dbl.BitsPerPixel = 64;
ref_u8           = ref; ref_u8.Datatype  = 'uint8';  ref_u8.BitsPerPixel  = 8;

niftiwrite(double(QSM),                  fullfile(output_dir, 'QSM.nii.gz'),         ref_dbl, 'Compressed', true);
niftiwrite(double(R2star),               fullfile(output_dir, 'R2star.nii.gz'),       ref_dbl, 'Compressed', true);
niftiwrite(uint8(Brain_Mask_vsharp),     fullfile(output_dir, 'brain_mask.nii.gz'),   ref_u8,  'Compressed', true);

mask_qs = logical(Brain_Mask_vsharp);
fprintf('QSM range:   [%.4f, %.4f] ppm\n', min(QSM(mask_qs)), max(QSM(mask_qs)));
fprintf('R2* range:   [%.1f, %.1f] s^-1\n', min(R2star(Brain_Mask)), max(R2star(Brain_Mask)));
disp('QSM reconstruction complete.')
fprintf('Output written to: %s\n', output_dir)
