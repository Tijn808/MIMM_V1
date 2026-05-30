% Sanity checks for MIMM pipeline outputs
% Runs 6 checks and saves a summary CSV + prints results to console.
%
% Checks:
%   1. WM-only stats (FA > 0.3 mask)
%   2. MVF vs FA atlas correlation (WM mask)
%   3. chi_myelin vs MVF correlation (WM mask)
%   4. Iron-rich region stats (R2* > 40 proxy mask)
%   5. Error map assessment
%   6. Basic vs Atlas MVF difference

paths_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'paths.m');
if ~exist(paths_file, 'file')
    error('paths.m not found. Copy MUMC_pipeline/paths_template.m to paths.m and fill in your paths.');
end
run(paths_file);
subj_dir = input_dir;
out_file  = fullfile(subj_dir, 'mimm', 'sanity_check_results.csv');

%% --- Load maps ---

disp('Loading maps...')
brain_mask  = logical(niftiread(fullfile(subj_dir, 'qsm',   'brain_mask.nii.gz')));
FA          = double(niftiread(fullfile(subj_dir, 'atlas',  'FA_atlas.nii.gz')));
R2star      = double(niftiread(fullfile(subj_dir, 'qsm',   'R2star.nii.gz')));
MVF_basic   = double(niftiread(fullfile(subj_dir, 'mimm',  'MVF_basic.nii.gz')));
MVF_atlas   = double(niftiread(fullfile(subj_dir, 'mimm',  'MVF_Atlas.nii.gz')));
FVF_basic   = double(niftiread(fullfile(subj_dir, 'mimm',  'FVF_basic.nii.gz')));
g_ratio     = double(niftiread(fullfile(subj_dir, 'mimm',  'g_ratio_basic.nii.gz')));
chi_myelin  = double(niftiread(fullfile(subj_dir, 'mimm',  'chi_myelin_basic.nii.gz')));
chi_iron    = double(niftiread(fullfile(subj_dir, 'mimm',  'chi_iron_est_basic.nii.gz')));
error_map   = double(niftiread(fullfile(subj_dir, 'mimm',  'error_basic.nii.gz')));

%% --- Masks ---

wm_mask   = brain_mask & (FA > 0.3);
iron_mask = brain_mask & (R2star > 40);
gm_mask   = brain_mask & (FA <= 0.3) & (FA > 0.05);

fprintf('\nBrain mask:  %d voxels\n', sum(brain_mask(:)));
fprintf('WM mask:     %d voxels (FA > 0.3)\n', sum(wm_mask(:)));
fprintf('GM mask:     %d voxels (0.05 < FA <= 0.3)\n', sum(gm_mask(:)));
fprintf('Iron mask:   %d voxels (R2* > 40)\n', sum(iron_mask(:)));

%% =========================================================
%  CHECK 1: WM-only statistics
% =========================================================

fprintf('\n========== CHECK 1: WM-ONLY STATS (FA > 0.3) ==========\n')

maps     = {MVF_basic, FVF_basic, g_ratio, R2star, chi_myelin, chi_iron};
mapnames = {'MVF', 'FVF', 'g_ratio', 'R2star', 'chi_myelin', 'chi_iron_est'};
units    = {'fraction', 'fraction', '-', 's-1', 'ppm', 'ppm'};

wm_stats = zeros(numel(maps), 4);  % mean, sd, p5, p95
for i = 1:numel(maps)
    vals = maps{i}(wm_mask);
    wm_stats(i,:) = [mean(vals), std(vals), prctile(vals,5), prctile(vals,95)];
    fprintf('%-15s  mean=%.4f  SD=%.4f  [P5=%.4f, P95=%.4f]  %s\n', ...
        mapnames{i}, wm_stats(i,1), wm_stats(i,2), wm_stats(i,3), wm_stats(i,4), units{i});
end

%% =========================================================
%  CHECK 2: MVF vs FA correlation in WM
% =========================================================

fprintf('\n========== CHECK 2: MVF vs FA CORRELATION (WM) ==========\n')

mvf_wm = MVF_basic(wm_mask);
fa_wm  = FA(wm_mask);
r_mvf_fa = pearson_r(mvf_wm, fa_wm);
fprintf('Pearson r(MVF, FA_atlas) in WM = %.4f\n', r_mvf_fa);
fprintf('Expected: positive correlation (both track WM integrity)\n');

%% =========================================================
%  CHECK 3: chi_myelin vs MVF correlation in WM
% =========================================================

fprintf('\n========== CHECK 3: chi_myelin vs MVF CORRELATION (WM) ==========\n')

chi_m_wm = chi_myelin(wm_mask);
r_mvf_chi = pearson_r(mvf_wm, chi_m_wm);
fprintf('Pearson r(MVF, chi_myelin) in WM = %.4f\n', r_mvf_chi);
fprintf('Expected: negative correlation (more myelin = more diamagnetic = more negative chi)\n');

%% =========================================================
%  CHECK 4: Iron-rich region stats
% =========================================================

fprintf('\n========== CHECK 4: IRON-RICH REGIONS (R2* > 40 s-1) ==========\n')

fprintf('R2*       iron regions:  mean=%.1f  SD=%.1f  [P5=%.1f, P95=%.1f] s-1\n', ...
    mean(R2star(iron_mask)), std(R2star(iron_mask)), ...
    prctile(R2star(iron_mask),5), prctile(R2star(iron_mask),95));
fprintf('chi_iron  iron regions:  mean=%.4f  SD=%.4f  [P5=%.4f, P95=%.4f] ppm\n', ...
    mean(chi_iron(iron_mask)), std(chi_iron(iron_mask)), ...
    prctile(chi_iron(iron_mask),5), prctile(chi_iron(iron_mask),95));
fprintf('chi_myelin iron regions: mean=%.4f  SD=%.4f  [P5=%.4f, P95=%.4f] ppm\n', ...
    mean(chi_myelin(iron_mask)), std(chi_myelin(iron_mask)), ...
    prctile(chi_myelin(iron_mask),5), prctile(chi_myelin(iron_mask),95));
fprintf('MVF       iron regions:  mean=%.4f  SD=%.4f\n', ...
    mean(MVF_basic(iron_mask)), std(MVF_basic(iron_mask)));
fprintf('Expected: chi_iron elevated, chi_myelin near 0 or less negative than WM\n');

% Compare iron-rich vs WM
fprintf('\nchi_iron  WM (FA>0.3):   mean=%.4f\n', mean(chi_iron(wm_mask)));
fprintf('chi_iron  iron regions:  mean=%.4f  (should be higher)\n', mean(chi_iron(iron_mask)));

%% =========================================================
%  CHECK 5: Error map
% =========================================================

fprintf('\n========== CHECK 5: ERROR MAP ==========\n')

fprintf('Error whole brain:  mean=%.4f  SD=%.4f  [P5=%.4f, P95=%.4f]\n', ...
    mean(error_map(brain_mask)), std(error_map(brain_mask)), ...
    prctile(error_map(brain_mask),5), prctile(error_map(brain_mask),95));
fprintf('Error WM (FA>0.3):  mean=%.4f  SD=%.4f\n', ...
    mean(error_map(wm_mask)), std(error_map(wm_mask)));
fprintf('Error GM (FA<0.3):  mean=%.4f  SD=%.4f\n', ...
    mean(error_map(gm_mask)), std(error_map(gm_mask)));
fprintf('Expected: lower error in WM (dictionary fits WM best)\n');

%% =========================================================
%  CHECK 6: Basic vs Atlas MVF difference
% =========================================================

fprintf('\n========== CHECK 6: BASIC vs ATLAS MVF ==========\n')

diff_map = MVF_basic - MVF_atlas;
fprintf('MVF difference (Basic - Atlas) whole brain:\n');
fprintf('  mean=%.4f  SD=%.4f  [P5=%.4f, P95=%.4f]\n', ...
    mean(diff_map(brain_mask)), std(diff_map(brain_mask)), ...
    prctile(diff_map(brain_mask),5), prctile(diff_map(brain_mask),95));
fprintf('MVF difference (Basic - Atlas) WM only:\n');
fprintf('  mean=%.4f  SD=%.4f  [P5=%.4f, P95=%.4f]\n', ...
    mean(diff_map(wm_mask)), std(diff_map(wm_mask)), ...
    prctile(diff_map(wm_mask),5), prctile(diff_map(wm_mask),95));
fprintf('Expected: small difference (atlas refines but should not dramatically change MVF)\n');
pct_changed = 100 * mean(abs(diff_map(wm_mask)) > 0.01);
fprintf('Voxels with |difference| > 0.01 in WM: %.1f%%\n', pct_changed);

%% =========================================================
%  Save CSV summary
% =========================================================

fid = fopen(out_file, 'w');
fprintf(fid, 'Check,Map,Mask,Mean,SD,P5,P95,Unit\n');

for i = 1:numel(maps)
    fprintf(fid, '1_WM_stats,%s,FA>0.3,%.4f,%.4f,%.4f,%.4f,%s\n', ...
        mapnames{i}, wm_stats(i,1), wm_stats(i,2), wm_stats(i,3), wm_stats(i,4), units{i});
end

fprintf(fid, '2_Correlation,r(MVF-FA),WM,%.4f,,,, \n', r_mvf_fa);
fprintf(fid, '3_Correlation,r(MVF-chiMyelin),WM,%.4f,,,, \n', r_mvf_chi);

fprintf(fid, '4_Iron,R2star,R2*>40,%.2f,%.2f,%.2f,%.2f,s-1\n', ...
    mean(R2star(iron_mask)), std(R2star(iron_mask)), prctile(R2star(iron_mask),5), prctile(R2star(iron_mask),95));
fprintf(fid, '4_Iron,chi_iron,R2*>40,%.4f,%.4f,%.4f,%.4f,ppm\n', ...
    mean(chi_iron(iron_mask)), std(chi_iron(iron_mask)), prctile(chi_iron(iron_mask),5), prctile(chi_iron(iron_mask),95));

fprintf(fid, '5_Error,error,brain,%.4f,%.4f,%.4f,%.4f,-\n', ...
    mean(error_map(brain_mask)), std(error_map(brain_mask)), prctile(error_map(brain_mask),5), prctile(error_map(brain_mask),95));
fprintf(fid, '5_Error,error,WM,%.4f,%.4f,%.4f,%.4f,-\n', ...
    mean(error_map(wm_mask)), std(error_map(wm_mask)), prctile(error_map(wm_mask),5), prctile(error_map(wm_mask),95));
fprintf(fid, '5_Error,error,GM,%.4f,%.4f,%.4f,%.4f,-\n', ...
    mean(error_map(gm_mask)), std(error_map(gm_mask)), prctile(error_map(gm_mask),5), prctile(error_map(gm_mask),95));

fprintf(fid, '6_BasicVsAtlas,MVF_diff,brain,%.4f,%.4f,%.4f,%.4f,fraction\n', ...
    mean(diff_map(brain_mask)), std(diff_map(brain_mask)), prctile(diff_map(brain_mask),5), prctile(diff_map(brain_mask),95));
fprintf(fid, '6_BasicVsAtlas,MVF_diff,WM,%.4f,%.4f,%.4f,%.4f,fraction\n', ...
    mean(diff_map(wm_mask)), std(diff_map(wm_mask)), prctile(diff_map(wm_mask),5), prctile(diff_map(wm_mask),95));

fclose(fid);
fprintf('\nResults saved to: %s\n', out_file);
disp('Sanity checks complete.')

function r = pearson_r(x, y)
    x = x(:); y = y(:);
    xm = x - mean(x); ym = y - mean(y);
    r = sum(xm .* ym) / sqrt(sum(xm.^2) * sum(ym.^2));
end
