% run_analysis.m — top-level analysis script for MUMC IMPROMYMS data
%
% Requires:
%   results_file : path to .mat with MIMM_basic, MIMM_DTI, MIMM_Atlas structs
%   mwf_file     : path to .mat with per-subject MWF values per ROI
%   chisep_file  : path to .mat with per-subject chi_neg values per ROI
%   clinical_file: path to .mat or .csv with clinical scores struct
%
% ROIs to analyse (update to match atlas labels used):
%   'NAWM', 'WM_lesion', 'whole_brain_WM'

%% --- Paths (adapt per study) ---

results_file  = '';   % e.g. 'Example_Results/stochastic_MIMM_results.mat'
mwf_file      = '';
chisep_file   = '';
clinical_file = '';
output_dir    = 'analysis_output';

if ~exist(output_dir, 'dir'); mkdir(output_dir); end

%% --- Load data ---

load(results_file,  'MIMM_basic', 'MIMM_DTI', 'MIMM_Atlas');
load(mwf_file,      'mwf_per_roi');       % struct: mwf_per_roi.NAWM, etc.
load(chisep_file,   'chisep_per_roi');    % struct: chisep_per_roi.NAWM, etc.
load(clinical_file, 'clinical');          % struct: clinical.edss, .t25fw, .nhpt, .ssdmt

%% --- Choose MIMM variant and ROIs ---

mimm_variants = {MIMM_basic, MIMM_DTI, MIMM_Atlas};
variant_names = {'Basic', 'DTI', 'Atlas'};
rois          = {'NAWM', 'WM_lesion'};

%% --- Run Bland-Altman and clinical correlation per variant and ROI ---

for vi = 1:numel(mimm_variants)
    mv = mimm_variants{vi};
    vname = variant_names{vi};

    for ri = 1:numel(rois)
        roi = rois{ri};
        out_sub = fullfile(output_dir, vname, roi);
        if ~exist(out_sub, 'dir'); mkdir(out_sub); end

        % Subject-level means within ROI — one value per subject.
        % These must be [N_subjects x 1] vectors, NOT voxel-level arrays.
        % Compute mean over ROI voxels for each subject before passing here.
        mvf_subj   = mv.MVF_roi.(roi);       % [N_subjects x 1] mean MVF per subject
        mwf_subj   = mwf_per_roi.(roi);      % [N_subjects x 1] mean MWF per subject
        chi_subj   = chisep_per_roi.(roi);   % [N_subjects x 1] mean |chi_neg| per subject

        %% Bland-Altman (subject-level: N points = N subjects)
        bland_altman(mvf_subj, mwf_subj, chi_subj, [vname '_' roi], out_sub);

        %% Clinical correlation
        params.mvf      = mv.MVF_roi.(roi);
        params.fvf      = mv.FVF_roi.(roi);
        params.g_ratio  = mv.g_ratio_roi.(roi);
        params.chi_iron = mv.chi_iron_roi.(roi);

        correlate_clinical(params, clinical, [vname '_' roi], out_sub);
    end
end

disp('Analysis complete.')
