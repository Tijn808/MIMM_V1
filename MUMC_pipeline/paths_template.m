% =========================================================================
% MUMC Pipeline — Machine Configuration Template
% =========================================================================
% 1. Copy this file to paths.m (same directory)
% 2. Fill in the four paths below
% 3. paths.m is gitignored — never committed
% =========================================================================

mimm_root  = '/path/to/MIMM';                % repo root (where RUNME.m lives)
input_dir  = '/path/to/raw/ME_GRE_data';     % folder with 501-ME_GRE_e*.nii.gz
output_dir = '/path/to/output/subject_id';   % created automatically if missing
chisep_dir = '/path/to/Chisep_Toolbox_v1.2.1';

% Derived subdirs — do not edit
qsm_dir    = fullfile(output_dir, 'qsm');
mimm_dir   = fullfile(output_dir, 'mimm');
chisep_out = fullfile(output_dir, 'chisep');
atlas_dir  = fullfile(output_dir, 'atlas');
fig_dir    = fullfile(output_dir, 'figures');
analysis_d = fullfile(output_dir, 'analysis');
sep_dir    = fullfile(output_dir, 'source_separation');
grase_dir  = fullfile(output_dir, 'grase');   % T2-GRASE MWF (registered to ME-GRE)
dti_dir    = fullfile(output_dir, 'dti');

for d = {qsm_dir, mimm_dir, chisep_out, atlas_dir, fig_dir, analysis_d, sep_dir, grase_dir, dti_dir}
    if ~exist(d{1}, 'dir'); mkdir(d{1}); end
end
