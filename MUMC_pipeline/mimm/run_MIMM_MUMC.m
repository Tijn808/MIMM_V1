% MIMM reconstruction for MUMC ME-GRE data
%
% Acquisition parameters (fixed for MUMC protocol):
%   Sequence : 3D ME-GRE (monopolar readout)
%   Scanner  : 3T Philips MR 7700
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
lambda_chi = 0.015;                      % QSM/magnitude weighting (paper L-curve default)

%% --- Paths (loaded from MUMC_pipeline/paths.m) ---

if ~exist('mimm_root', 'var')
    paths_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'paths.m');
    if ~exist(paths_file, 'file')
        error('paths.m not found. Copy MUMC_pipeline/paths_template.m to paths.m and fill in your paths.');
    end
    run(paths_file);
end
output_dir = mimm_dir;

% lambda_chi precedence (each lives in a data file, not hard-coded, so the
% choice stays explicit and reversible — delete the file to fall back):
%   1. cohort-locked value (next to paths.m)   — shared by all subjects
%   2. per-subject sweep recommendation        — this subject's own knee
%   3. paper default 0.015                      — if neither file exists
cohort_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'lambda_chi_cohort.txt');
rec_file    = fullfile(mimm_dir, 'lambda_chi_recommended.txt');
if exist(cohort_file, 'file')
    lambda_chi = str2double(strtrim(fileread(cohort_file)));
    fprintf('Using cohort-locked lambda_chi = %.4g (from %s).\n', lambda_chi, cohort_file);
elseif exist(rec_file, 'file')
    lambda_chi = str2double(strtrim(fileread(rec_file)));
    fprintf('Using per-subject lambda_chi = %.4g (from %s).\n', lambda_chi, rec_file);
else
    fprintf('Using default lambda_chi = %.4g (run sweep_lambda_chi.m to tune).\n', lambda_chi);
end

%% --- Add MIMM toolbox to path ---

run(fullfile(mimm_root, 'MIMM_set_path.m'));

%% --- Load dictionary ---

% Use MUMC-specific dictionary if available (no interpolation needed).
% Fall back to original paper dictionary with interpolation if not yet generated.
mumc_dict = fullfile(mimm_root, 'Dictionary', 'MIMM_dictionary_stochastic_MUMC.mat');
orig_dict = fullfile(mimm_root, 'Dictionary', 'MIMM_dictionary_stochastic.mat');
if exist(mumc_dict, 'file')
    dict_file = mumc_dict;
    disp('Using MUMC dictionary (TEs: 6/12/18/24/30 ms — no interpolation).');
else
    dict_file = orig_dict;
    warning('MUMC dictionary not found — using original with TE interpolation.');
end
stoch = load(dict_file, 'dictionary');
dict  = stoch.dictionary;

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
theta_file_atlas = fullfile(subj_dir, 'atlas', 'theta_atlas.nii.gz');   % degrees, from register_atlas.sh

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
