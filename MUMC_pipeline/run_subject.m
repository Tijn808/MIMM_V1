function run_subject(subj_dir_arg, mimm_root_arg, chisep_dir_arg, steps_arg, prefix_arg)
% Run the MATLAB pipeline steps for one subject without needing paths.m.
%
% Called by run_cohort.sh for each subject. Sets all workspace variables
% that paths.m would set, then runs each pipeline script. Scripts detect
% the pre-set variables and skip their own paths.m loading (mimm_root guard).
%
% Usage:
%   run_subject(subj_dir, mimm_root, chisep_dir)
%   run_subject(subj_dir, mimm_root, chisep_dir, {'prepare_mgre','qsm'})
%   run_subject(subj_dir, mimm_root, chisep_dir, {'prepare_mgre','qsm'}, '501')
%
% steps_arg (optional cell array):
%   'prepare_mgre'  — prepare_mgre.m  (default first half, run before register_atlas)
%   'qsm'           — run_QSM_chisep.m  (writes mag_e1.nii.gz needed by register)
%   'chisep'        — run_chisep_MUMC.m (default second half, run after register)
%   'mimm'          — run_MIMM_MUMC.m
%   'grase'         — grase/ingest_grase.m (skipped if MWF_native.nii.gz absent)
%
% prefix_arg (optional string):
%   Series number prefix for prepare_mgre.m (e.g. '501').
%   Auto-detected from the first *-ME_GRE_e1.nii.gz file if not supplied.
%
% IMPORTANT — correct call order from run_cohort.sh:
%   1. run_subject(..., {'prepare_mgre','qsm'})   % creates mag_e1.nii.gz
%   2. register_atlas.sh                          % needs mag_e1.nii.gz
%   3. preprocess_dti.sh                          % needs mag_e1.nii.gz
%   4. run_subject(..., {'chisep','mimm'})         % needs atlas outputs

if nargin < 4 || isempty(steps_arg)
    steps_arg = {'prepare_mgre', 'qsm', 'chisep', 'mimm'};
end
if nargin < 5; prefix_arg = ''; end

% ── Set all path variables (mirrors paths_template.m) ────────────────────────
% Variables set here are consumed by the scripts called via run() below.
% MATLAB's analyser cannot see into run() callers so NASGU warnings are expected.
input_dir  = subj_dir_arg;
output_dir = subj_dir_arg;
subj_dir   = subj_dir_arg;   %#ok<NASGU> % alias used by run_MIMM_MUMC.m and others
mimm_root  = mimm_root_arg;  % triggers the paths.m guard in each pipeline script
chisep_dir = chisep_dir_arg; %#ok<NASGU>

qsm_dir    = fullfile(output_dir, 'qsm');
mimm_dir   = fullfile(output_dir, 'mimm');
chisep_out = fullfile(output_dir, 'chisep');
atlas_dir  = fullfile(output_dir, 'atlas');
fig_dir    = fullfile(output_dir, 'figures');
analysis_d = fullfile(output_dir, 'analysis');
sep_dir    = fullfile(output_dir, 'source_separation');
grase_dir  = fullfile(output_dir, 'grase');
dti_dir    = fullfile(output_dir, 'dti');

for d = {qsm_dir, mimm_dir, chisep_out, atlas_dir, fig_dir, analysis_d, ...
         sep_dir, grase_dir, dti_dir}
    if ~exist(d{1}, 'dir'); mkdir(d{1}); end
end

% Read cohort-level lambda_chi if locked (see sweep_lambda_chi.m)
lambda_chi = 0.015;
cohort_lam = fullfile(fileparts(mfilename('fullpath')), 'lambda_chi_cohort.txt');
rec_lam    = fullfile(mimm_dir, 'lambda_chi_recommended.txt');
if exist(cohort_lam, 'file')
    lambda_chi = str2double(strtrim(fileread(cohort_lam)));
    fprintf('[run_subject] lambda_chi = %.4g (cohort-locked)\n', lambda_chi);
elseif exist(rec_lam, 'file')
    lambda_chi = str2double(strtrim(fileread(rec_lam)));
    fprintf('[run_subject] lambda_chi = %.4g (per-subject sweep)\n', lambda_chi);
else
    fprintf('[run_subject] lambda_chi = %.4g (default)\n', lambda_chi);
end

% Add MIMM toolbox to path
run(fullfile(mimm_root, 'MIMM_set_path.m'));

pipeline_dir = fileparts(mfilename('fullpath'));

% ── Run requested steps ───────────────────────────────────────────────────────
for s = steps_arg
    step = s{1};
    fprintf('\n[run_subject] %s — step: %s\n', subj_dir_arg, step);
    t0 = tic;

    switch step
        case 'prepare_mgre'
            % Resolve the series-number prefix. Use supplied value if given;
            % otherwise auto-detect from the first *-ME_GRE_e1.nii.gz in input_dir.
            has_mumc = exist(fullfile(input_dir, 'gremag.nii'),    'file') || ...
                       exist(fullfile(input_dir, 'gremag.nii.gz'), 'file');
            if ~isempty(prefix_arg)
                prefix = prefix_arg; %#ok<NASGU>
            elseif has_mumc
                prefix = ''; %#ok<NASGU>  % MUMC 4D format; prepare_mgre auto-detects it
                fprintf('  [MUMC format] gremag/grepha detected\n');
            else
                hits = dir(fullfile(input_dir, '*-ME_GRE_e1.nii.gz'));
                if isempty(hits)
                    error('[run_subject] No *-ME_GRE_e1.nii.gz found in %s. Supply prefix_arg.', input_dir);
                end
                prefix = strrep(hits(1).name, '-ME_GRE_e1.nii.gz', '');
                fprintf('  [auto-prefix] detected prefix: %s\n', prefix);
            end
            run(fullfile(pipeline_dir, 'preprocess', 'prepare_mgre.m'));

        case 'qsm'
            run(fullfile(pipeline_dir, 'qsm', 'run_QSM_chisep.m'));

        case 'chisep'
            run(fullfile(pipeline_dir, 'qsm', 'run_chisep_MUMC.m'));

        case 'mimm'
            run(fullfile(pipeline_dir, 'mimm', 'run_MIMM_MUMC.m'));

        case 'grase'
            % Run if any accepted native MWF map is present (MUMC outputs MWF.nii).
            have_mwf = exist(fullfile(grase_dir, 'MWF_native.nii.gz'), 'file') || ...
                       exist(fullfile(grase_dir, 'MWF_native.nii'),    'file') || ...
                       exist(fullfile(grase_dir, 'MWF.nii'),           'file');
            if have_mwf
                run(fullfile(pipeline_dir, 'grase', 'ingest_grase.m'));
            else
                fprintf('  [skip] no native MWF map in %s (MWF_native.nii.gz / MWF.nii)\n', grase_dir);
            end

        otherwise
            warning('[run_subject] Unknown step: %s — skipped.', step);
    end

    fprintf('[run_subject] %s done in %.1f min\n', step, toc(t0)/60);
end

fprintf('\n[run_subject] All steps complete for: %s\n', subj_dir_arg);
end
