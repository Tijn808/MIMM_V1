% Lesion segmentation using LST-LPA (Lesion Prediction Algorithm).
%
% LST-LPA uses T1w + T2-FLAIR to segment WM lesions. It requires SPM12
% and the LST toolbox (https://www.applied-statistics.de/lst.html).
%
% Inputs expected in <subj_dir>/lesion/:
%   FLAIR_native.nii.gz    T2-FLAIR in native space (uncompressed copy created)
%
% Inputs expected in <subj_dir>/t1w/:
%   T1w_native.nii.gz      T1-weighted image in native space
%
% Output:
%   lesion/lesion_mask_native.nii.gz   Binary lesion mask in FLAIR space
%
% After this script, run register_flair.sh to bring the mask to ME-GRE space.
%
% Usage:
%   Fill in paths.m, then run this script from MATLAB.
%   LST must be installed: spm_get_defaults('LST') should not error.

%% --- Paths ---
if ~exist('mimm_root', 'var')
    paths_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'paths.m');
    if ~exist(paths_file, 'file')
        error('paths.m not found. Copy MUMC_pipeline/paths_template.m to paths.m.');
    end
    run(paths_file);
end

lesion_d = fullfile(output_dir, 'lesion');
t1w_d    = fullfile(output_dir, 't1w');

flair_gz = fullfile(lesion_d, 'FLAIR_native.nii.gz');
t1w_gz   = fullfile(t1w_d,    'T1w_native.nii.gz');

if ~exist(flair_gz, 'file')
    error('FLAIR_native.nii.gz not found in %s.\nCopy the T2-FLAIR NIfTI there first.', lesion_d);
end

%% --- Check SPM + LST are available (auto-add from known location) ---
% Honour $SPM_DIR if set, else the local install at ~/spm12.
if ~exist('spm', 'file')
    spm_dir = getenv('SPM_DIR');
    if isempty(spm_dir); spm_dir = fullfile(getenv('HOME'), 'spm12'); end
    if exist(fullfile(spm_dir, 'spm.m'), 'file')
        addpath(spm_dir);
        addpath(fullfile(spm_dir, 'toolbox', 'LST'));
    else
        error(['SPM12 not found. Install it (e.g. git clone ' ...
               'https://github.com/spm/spm12 ~/spm12) or set $SPM_DIR, ' ...
               'then re-run.']);
    end
end
try
    spm('Defaults', 'fMRI');
catch
    error('SPM12 initialisation failed.');
end
if ~exist('ps_LST_lpa', 'file')
    % SPM is on path but LST toolbox isn't — try the standard toolbox location.
    lst_dir = fullfile(fileparts(which('spm')), 'toolbox', 'LST');
    if exist(fullfile(lst_dir, 'ps_LST_lpa.m'), 'file')
        addpath(lst_dir);
    else
        error(['LST toolbox not found. Install LST_3.0.0.zip from ' ...
               'https://www.applied-statistics.de/lst.html into ' ...
               'spm12/toolbox/LST.']);
    end
end

%% --- Decompress inputs (LST requires uncompressed .nii) ---
fprintf('Decompressing inputs...\n');
flair_nii = gunzip(flair_gz, lesion_d); flair_nii = flair_nii{1};

if exist(t1w_gz, 'file')
    t1w_nii = gunzip(t1w_gz, t1w_d); t1w_nii = t1w_nii{1};
    use_t1w = true;
    fprintf('T1w found — running LST-LPA with T1w + FLAIR (recommended).\n');
else
    t1w_nii = '';
    use_t1w = false;
    fprintf('T1w not found — running LST-LPA with FLAIR only.\n');
end

%% --- Run LST-LPA ---
fprintf('Running LST-LPA lesion segmentation...\n');
if use_t1w
    ps_LST_lpa(flair_nii, t1w_nii);
else
    ps_LST_lpa(flair_nii);
end

%% --- Locate and threshold output ---
% LST-LPA writes a probability map named ples_lpa_<prefix>.nii
[lesion_dir_str, flair_base] = fileparts(flair_nii);
ples_file = fullfile(lesion_dir_str, ['ples_lpa_' flair_base '.nii']);
if ~exist(ples_file, 'file')
    % LST may prefix with 'm' when T1w was provided
    ples_file = fullfile(lesion_dir_str, ['ples_lpa_m' flair_base '.nii']);
end
if ~exist(ples_file, 'file')
    error('LST output not found. Expected: %s', ples_file);
end

% Threshold at p > 0.5 to get binary mask
ples_vol = double(niftiread(ples_file));
mask_vol = uint8(ples_vol > 0.5);

% Save as NIfTI
ref      = niftiinfo(ples_file);
ref.Datatype    = 'uint8';
ref.BitsPerPixel = 8;
mask_out = fullfile(lesion_d, 'lesion_mask_native.nii.gz');
niftiwrite(mask_vol, mask_out, ref, 'Compressed', true);

n_lesion = sum(mask_vol(:));
fprintf('Lesion mask written: %s\n', mask_out);
fprintf('  Lesion voxels (p>0.5): %d  (%.1f mL at 1mm iso)\n', n_lesion, n_lesion/1000);

%% --- Clean up uncompressed intermediates ---
delete(flair_nii);
if use_t1w && exist(t1w_nii, 'file'); delete(t1w_nii); end
delete(ples_file);

fprintf('\nNext step: bash register_flair.sh %s\n', output_dir);
