function run_subject(subj_dir_arg, mimm_root_arg, chisep_dir_arg, steps_arg)
% Run the MATLAB pipeline steps for one subject without needing paths.m.
%
% Called by run_cohort.sh for each subject. Sets all workspace variables
% that paths.m would set, then runs each pipeline script. Scripts detect
% the pre-set variables and skip their own paths.m loading (mimm_root guard).
%
% Usage (from MATLAB command line or -batch):
%   run_subject('/path/to/cohort/sub-01', '/path/to/MIMM', '/path/to/chisep')
%   run_subject('/path/to/sub-01', '/path/to/MIMM', '/path/to/chisep', ...
%               {'prepare_mgre','qsm','chisep','mimm'})
%
% steps_arg (optional cell array, default = all):
%   'prepare_mgre'  — prepare_mgre.m   (requires raw NIfTIs in subj_dir)
%   'qsm'           — run_QSM_chisep.m
%   'chisep'        — run_chisep_MUMC.m
%   'mimm'          — run_MIMM_MUMC.m
%   'grase'         — grase/ingest_grase.m  (skipped if MWF_native not present)
%
% Note: shell-based steps (register_atlas.sh, preprocess_dti.sh) are run
% separately by run_cohort.sh before this function is called.

if nargin < 4 || isempty(steps_arg)
    steps_arg = {'prepare_mgre', 'qsm', 'chisep', 'mimm'};
end

% ── Set all path variables (mirrors paths_template.m) ────────────────────────
input_dir  = subj_dir_arg;
output_dir = subj_dir_arg;
mimm_root  = mimm_root_arg;   % this variable being set triggers the guard in each script
chisep_dir = chisep_dir_arg;

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
            run(fullfile(pipeline_dir, 'preprocess', 'prepare_mgre.m'));

        case 'qsm'
            run(fullfile(pipeline_dir, 'qsm', 'run_QSM_chisep.m'));

        case 'chisep'
            run(fullfile(pipeline_dir, 'qsm', 'run_chisep_MUMC.m'));

        case 'mimm'
            run(fullfile(pipeline_dir, 'mimm', 'run_MIMM_MUMC.m'));

        case 'grase'
            mwf_native = fullfile(grase_dir, 'MWF_native.nii.gz');
            if exist(mwf_native, 'file')
                run(fullfile(pipeline_dir, 'grase', 'ingest_grase.m'));
            else
                fprintf('  [skip] MWF_native.nii.gz not found in %s\n', grase_dir);
            end

        otherwise
            warning('[run_subject] Unknown step: %s — skipped.', step);
    end

    fprintf('[run_subject] %s done in %.1f min\n', step, toc(t0)/60);
end

fprintf('\n[run_subject] All steps complete for: %s\n', subj_dir_arg);
end
