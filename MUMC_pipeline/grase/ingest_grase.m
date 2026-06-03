% Register T2-GRASE MWF map to ME-GRE space for each subject.
%
% The T2-GRASE scan is acquired at 1.5x1.5x4 mm (different from ME-GRE 1 mm iso),
% so it needs to be registered and resliced before it can be compared voxel-for-
% voxel with MIMM MVF maps. Registration uses a rigid-body (6 DOF) flirt to align
% the GRASE volume to the ME-GRE echo 1 magnitude.
%
% Inputs expected in subj_dir/grase/  (native MWF map — first match used):
%   MWF_native.nii.gz / MWF_native.nii / MWF.nii   MWF map in T2-GRASE space.
%   MUMC's MWFfit.m (020_MWF_calc.sh) outputs an uncompressed MWF.nii — drop it
%   in directly, no rename needed (the pipeline OUTPUT is MWF.nii.gz, distinct).
%
% Registration source (best contrast first — first match used):
%   b0.nii.gz / b0.nii                 explicit registration target, OR
%   grase.nii.gz / grase.nii           4D GRASE series → first volume extracted, OR
%   (fallback) the MWF map itself      low contrast, last resort
%
% Outputs (written to subj_dir/grase/):
%   MWF.nii.gz          — MWF in ME-GRE space (1 mm iso, same grid as MVF maps)
%
% Usage:
%   paths.m must be filled in. Run from MATLAB:
%     run('MUMC_pipeline/grase/ingest_grase.m')

%% --- Paths ---
if ~exist('mimm_root', 'var')
    paths_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'paths.m');
    if ~exist(paths_file, 'file')
        error('paths.m not found. Copy MUMC_pipeline/paths_template.m to paths.m and fill in your paths.');
    end
    run(paths_file);
end

% FSL path (same logic as preprocess_dti.sh / register_atlas.sh)
if ~isempty(getenv('FSLDIR'))
    fsl = fullfile(getenv('FSLDIR'), 'bin');
else
    fsl = '/home/tijn-saes/fsl/bin';
end

out_dir = grase_dir;   % subj_dir/grase/

%% --- Locate the native MWF map (accept MUMC and pipeline naming) ---
ref = fullfile(qsm_dir, 'mag_e1.nii.gz');
if ~exist(ref, 'file')
    error('ME-GRE echo-1 reference not found: %s\nRun run_QSM_chisep.m first.', ref);
end

mwf_candidates = {'MWF_native.nii.gz', 'MWF_native.nii', 'MWF.nii'};
mwf_native = '';
for k = 1:numel(mwf_candidates)
    cand = fullfile(out_dir, mwf_candidates{k});
    if exist(cand, 'file'); mwf_native = cand; break; end
end
if isempty(mwf_native)
    if exist(fullfile(out_dir, 'MWF.nii.gz'), 'file')
        error(['Found MWF.nii.gz in %s, but that is the pipeline OUTPUT name.\n' ...
               'Rename the native MUMC map to MWF_native.nii.gz to avoid a collision.'], out_dir);
    end
    error(['Native MWF map not found in %s.\n' ...
           'Accepted: MWF_native.nii.gz, MWF_native.nii, or MWF.nii (MUMC output).'], out_dir);
end
fprintf('Native MWF map: %s\n', mwf_native);

%% --- Choose registration source (best contrast first) ---
reg_src    = '';
grase_vol0 = '';   % tracked for cleanup if we extract it

% 1. explicit b0 / registration target
for c = {'b0.nii.gz', 'b0.nii'}
    cand = fullfile(out_dir, c{1});
    if exist(cand, 'file'); reg_src = cand; fprintf('Registration source: %s\n', c{1}); break; end
end

% 2. first volume of the 4D GRASE series (T2-weighted — good registration target)
if isempty(reg_src)
    for c = {'grase.nii.gz', 'grase.nii'}
        g = fullfile(out_dir, c{1});
        if exist(g, 'file')
            grase_vol0 = fullfile(out_dir, 'grase_vol0.nii.gz');
            st = system(sprintf('"%s/fslroi" "%s" "%s" 0 1', fsl, g, grase_vol0));
            if st ~= 0; error('fslroi failed extracting first GRASE volume from %s.', c{1}); end
            reg_src = grase_vol0;
            fprintf('Registration source: first volume of %s\n', c{1});
            break;
        end
    end
end

% 3. fall back to the MWF map itself (low contrast, last resort)
if isempty(reg_src)
    reg_src = mwf_native;
    disp('No b0/grase found — registering the MWF map directly (lower contrast).');
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

%% --- Clean up temporary registration volume ---
if ~isempty(grase_vol0) && exist(grase_vol0, 'file')
    delete(grase_vol0);
end

%% --- Unit check: ensure MWF is a fraction (0-1), not a percentage (0-100) ---
% MIMM MVF is a fraction; the comparison figures assume MWF is too. MUMC's
% MWFfit may return percent. If the in-brain median exceeds 1.0, rescale by /100.
mwf_vol   = double(niftiread(mwf_out));
brain_vol = logical(niftiread(fullfile(qsm_dir, 'brain_mask.nii.gz')));
wm_vals   = mwf_vol(brain_vol & (mwf_vol > 0));

if ~isempty(wm_vals) && median(wm_vals) > 1.0
    fprintf('  MWF median = %.1f (>1) — looks like percent; rescaling /100 to fraction.\n', ...
        median(wm_vals));
    mwf_vol = mwf_vol / 100;
    out_info = niftiinfo(mwf_out);
    niftiwrite(mwf_vol, mwf_out, out_info, 'Compressed', true);
    wm_vals = mwf_vol(brain_vol & (mwf_vol > 0));
end

%% --- Sanity check ---
info = niftiinfo(mwf_out);
fprintf('\nDone. MWF.nii.gz written to: %s\n', mwf_out);
fprintf('  Grid: %s  VoxelSize: %s mm\n', mat2str(info.ImageSize), mat2str(info.PixelDimensions));

if ~isempty(wm_vals)
    fprintf('  MWF brain percentiles  P5=%.3f  P50=%.3f  P95=%.3f\n', ...
        prctile(wm_vals, 5), prctile(wm_vals, 50), prctile(wm_vals, 95));
    fprintf('  Expected WM range: ~0.05-0.25 (McCreary et al., Laule et al.).\n');
end
