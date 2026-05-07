% MIMM reconstruction for MUMC ME-GRE data
%
% Acquisition parameters (fixed for MUMC protocol):
%   Sequence : 3D ME-GRE (monopolar readout)
%   Scanner  : 3T Siemens
%   TE       : 6, 12, 18, 24, 30 ms
%   TR       : 35 ms
%   Voxel    : 1 x 1 x 1 mm isotropic
%   Matrix   : 224 x 224 x 154
%
% Prerequisites (must be completed before this script):
%   1. QSM reconstruction       — MUMC_pipeline/qsm/run_QSM.m
%   2. DTI preprocessing        — topup + eddy + dtifit (or use atlas theta)
%   3. All outputs registered to ME-GRE space
%
% Inputs expected in subj_dir:
%   magnitude.nii.gz        4D ME-GRE magnitude  [224 x 224 x 154 x 5]
%   qsm/QSM.nii.gz          QSM map              [224 x 224 x 154]  (ppm)
%   qsm/brain_mask.nii.gz   Binary brain mask    [224 x 224 x 154]
%   dti/FA.nii.gz           FA map               [224 x 224 x 154]  (optional)
%   dti/theta.nii.gz        Fibre angle map      [224 x 224 x 154]  (optional, degrees)
%   atlas/FA_atlas.nii.gz   Atlas FA             [224 x 224 x 154]  (optional)
%   atlas/theta_atlas.nii.gz Atlas theta         [224 x 224 x 154]  (optional, degrees)

%% --- Acquisition parameters (do not change for MUMC data) ---

TE     = [6, 12, 18, 24, 30] * 1e-3;   % echo times in seconds
lambda_chi = 0.015;                      % QSM/magnitude weighting (L-curve optimised)

%% --- Paths (adapt per subject) ---

mimm_root = '/home/tijn-saes/Documents/Internship/MIMM';
subj_dir  = '';      % e.g. '/data/MUMC/sub-01/'
output_dir = fullfile(subj_dir, 'mimm');

if ~exist(output_dir, 'dir'); mkdir(output_dir); end

%% --- Add MIMM toolbox to path ---

run(fullfile(mimm_root, 'MIMM_set_path.m'));

%% --- Load dictionary ---

dict_file = fullfile(mimm_root, 'Dictionary', 'MIMM_dictionary_stochastic.mat');
stoch     = load(dict_file, 'dictionary');
dict      = stoch.dictionary;

%% --- Load input data ---

mag_data   = double(niftiread(fullfile(subj_dir, 'magnitude.nii.gz')));   % [224x224x154x5]
QSM        = double(niftiread(fullfile(subj_dir, 'qsm', 'QSM.nii.gz'))); % [224x224x154]
Brain_Mask = logical(niftiread(fullfile(subj_dir, 'qsm', 'brain_mask.nii.gz')));

% MIMM expects iField = complex magnitude (real-valued ok; imaginary = 0)
% Normalise each echo relative to first echo to form decay curves
iField = mag_data;   % [224 x 224 x 154 x 5], real-valued, already positive

%% --- Optional: load orientation maps ---

dti_available   = false;
atlas_available = false;

fa_file_dti      = fullfile(subj_dir, 'dti',   'FA.nii.gz');
theta_file_dti   = fullfile(subj_dir, 'dti',   'theta.nii.gz');
fa_file_atlas    = fullfile(subj_dir, 'atlas', 'FA_atlas.nii.gz');
theta_file_atlas = fullfile(subj_dir, 'atlas', 'theta_atlas.nii.gz');

if exist(fa_file_dti, 'file') && exist(theta_file_dti, 'file')
    FA_DTI    = double(niftiread(fa_file_dti));
    theta_DTI = double(niftiread(theta_file_dti));
    dti_available = true;
    disp('DTI orientation maps loaded.');
end

if exist(fa_file_atlas, 'file') && exist(theta_file_atlas, 'file')
    FA_atlas    = double(niftiread(fa_file_atlas));
    theta_atlas = double(niftiread(theta_file_atlas));
    atlas_available = true;
    disp('Atlas orientation maps loaded.');
end

%% --- Run MIMM ---

% Basic MIMM (no orientation prior)
disp('Running MIMM Basic...');
MIMM_basic = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, 'basic');
disp('Basic done.');

% DTI orientation-informed MIMM
if dti_available
    disp('Running MIMM DTI...');
    MIMM_DTI = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, ...
                    'orientation_informed', FA_DTI, theta_DTI);
    disp('DTI done.');
else
    MIMM_DTI = [];
    warning('DTI maps not found — MIMM_DTI skipped.');
end

% Atlas orientation-informed MIMM
if atlas_available
    disp('Running MIMM Atlas...');
    MIMM_Atlas = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, ...
                      'orientation_informed', FA_atlas, theta_atlas);
    disp('Atlas done.');
else
    MIMM_Atlas = [];
    warning('Atlas maps not found — MIMM_Atlas skipped.');
end

%% --- Save results ---

save_file = fullfile(output_dir, 'MIMM_results.mat');
save(save_file, 'MIMM_basic', 'MIMM_DTI', 'MIMM_Atlas');
disp(['Results saved to: ' save_file]);

%% --- Export output maps as NIfTI ---
% Uses QSM header as spatial reference (same grid)

ref_info = niftiinfo(fullfile(subj_dir, 'qsm', 'QSM.nii.gz'));
export_nifti(MIMM_basic, 'basic',  ref_info, output_dir);
if dti_available;   export_nifti(MIMM_DTI,   'DTI',   ref_info, output_dir); end
if atlas_available; export_nifti(MIMM_Atlas,  'Atlas', ref_info, output_dir); end

disp('MIMM reconstruction complete.');


%% --- Helper: write MIMM output fields as NIfTI ---
function export_nifti(mimm, tag, ref_info, out_dir)
    if isempty(mimm); return; end
    fields = {'MVF', 'FVF', 'g_ratio', 'R2s', 'chi_iron_est', 'chi_myelin', ...
              'theta_est', 'error'};
    for fi = 1:numel(fields)
        fn   = fields{fi};
        data = mimm.(fn);
        info = ref_info;
        info.Datatype        = 'double';
        info.ImageSize       = size(data);
        info.PixelDimensions = ref_info.PixelDimensions(1:3);
        fname = fullfile(out_dir, sprintf('%s_%s.nii.gz', fn, tag));
        niftiwrite(data, fname, info, 'Compressed', true);
    end
end
