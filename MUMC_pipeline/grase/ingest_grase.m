% Register T2-GRASE MWF map to ME-GRE space for each subject.
%
% The T2-GRASE scan is acquired at 1.5x1.5x4 mm (different from ME-GRE 1 mm iso),
% so it needs to be registered and resliced before it can be compared voxel-for-
% voxel with MIMM MVF maps. Registration uses a rigid-body (6 DOF) flirt to align
% the GRASE volume to the ME-GRE echo 1 magnitude.
%
% Inputs expected in grase_raw_dir:
%   MWF_native.nii.gz   — MWF map in native T2-GRASE space (from T2-relaxometry)
%   b0.nii.gz           — any GRASE volume suitable for registration (b=0 or T2w)
%                         (if absent, MWF_native is used directly as the reference)
%
% Outputs (written to subj_dir/grase/):
%   MWF.nii.gz          — MWF in ME-GRE space (1 mm iso, same grid as MVF maps)
%
% Usage:
%   paths.m must be filled in. Run from MATLAB:
%     run('MUMC_pipeline/grase/ingest_grase.m')

%% --- Paths ---
paths_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'paths.m');
if ~exist(paths_file, 'file')
    error('paths.m not found. Copy MUMC_pipeline/paths_template.m to paths.m and fill in your paths.');
end
run(paths_file);

% FSL path (same logic as preprocess_dti.sh / register_atlas.sh)
if ~isempty(getenv('FSLDIR'))
    fsl = fullfile(getenv('FSLDIR'), 'bin');
else
    fsl = '/home/tijn-saes/fsl/bin';
end

out_dir = grase_dir;   % subj_dir/grase/

%% --- Locate inputs ---
mwf_native = fullfile(out_dir, 'MWF_native.nii.gz');
b0_grase   = fullfile(out_dir, 'b0.nii.gz');
ref        = fullfile(qsm_dir, 'mag_e1.nii.gz');

if ~exist(mwf_native, 'file')
    error('MWF_native.nii.gz not found in %s.\nCopy the T2-GRASE MWF map there first.', out_dir);
end
if ~exist(ref, 'file')
    error('ME-GRE echo-1 reference not found: %s\nRun run_QSM_chisep.m first.', ref);
end

% Use b0_grase as registration source if present; otherwise register MWF directly.
if exist(b0_grase, 'file')
    reg_src = b0_grase;
    disp('Using b0.nii.gz as registration source.');
else
    reg_src = mwf_native;
    disp('b0.nii.gz not found — registering MWF_native.nii.gz directly to ME-GRE.');
end

%% --- Rigid-body registration: GRASE → ME-GRE ---
omat = fullfile(out_dir, 'grase2mgre.mat');
reg_out = fullfile(out_dir, 'b0_in_mgre.nii.gz');

fprintf('Running flirt (GRASE → ME-GRE, 6 DOF)...\n');
cmd = sprintf('"%s/flirt" -in "%s" -ref "%s" -omat "%s" -out "%s" -dof 6 -cost normcorr', ...
    fsl, reg_src, ref, omat, reg_out);
status = system(cmd);
if status ~= 0; error('flirt registration failed.'); end

%% --- Apply transform to MWF map ---
mwf_out = fullfile(out_dir, 'MWF.nii.gz');
fprintf('Applying transform to MWF map...\n');
cmd = sprintf('"%s/flirt" -in "%s" -ref "%s" -applyxfm -init "%s" -out "%s" -interp trilinear', ...
    fsl, mwf_native, ref, omat, mwf_out);
status = system(cmd);
if status ~= 0; error('flirt applyxfm failed.'); end

%% --- Sanity check ---
info = niftiinfo(mwf_out);
fprintf('\nDone. MWF.nii.gz written to: %s\n', mwf_out);
fprintf('  Grid: %s  VoxelSize: %s mm\n', mat2str(info.ImageSize), mat2str(info.PixelDimensions));

mwf_vol   = double(niftiread(mwf_out));
brain_vol = logical(niftiread(fullfile(qsm_dir, 'brain_mask.nii.gz')));
wm_vals   = mwf_vol(brain_vol & (mwf_vol > 0));
if ~isempty(wm_vals)
    fprintf('  MWF brain percentiles  P5=%.3f  P50=%.3f  P95=%.3f\n', ...
        prctile(wm_vals, 5), prctile(wm_vals, 50), prctile(wm_vals, 95));
    fprintf('  Expected WM range: ~0.05-0.25 (McCreary et al., Laule et al.).\n');
end
