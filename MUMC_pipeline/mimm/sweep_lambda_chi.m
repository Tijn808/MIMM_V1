% Sweep lambda_chi to diagnose the MIMM chi_total vs QSM bias (figure 30).
%
% lambda_chi weights the QSM-matching term against the magnitude-decay term
% in MIMM's voxel-wise cost (MIMM.m line 123):
%     min over dict of:  lambda_chi*|chi_total_dict - QSM| + (1 - mag_corr)
% A larger lambda_chi forces chi_total closer to the measured QSM, at the cost
% of a looser magnitude-curve fit. This script finds how slope/bias/r of
% chi_total-vs-QSM respond to lambda_chi, so the value can be chosen from data
% rather than inherited from the paper's L-curve.
%
% Efficiency: the magnitude correlation (the expensive dictionary x data
% matmul) and chi_error are computed ONCE per slice; only the cheap argmin is
% repeated per lambda. A full sweep costs ~one MIMM run.
%
% Basic strategy only (lambda_theta = 0) to isolate the lambda_chi effect.
%
% Outputs (to mimm_dir):
%   lambda_chi_sweep.csv   columns: lambda_chi, bias, slope, r, mean_abs_err
%   lambda_chi_sweep.png   slope / bias / r vs lambda_chi, default marked

%% --- Paths ---
if ~exist('mimm_root', 'var')
    paths_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'paths.m');
    if ~exist(paths_file, 'file')
        error('paths.m not found. Copy MUMC_pipeline/paths_template.m to paths.m and fill in your paths.');
    end
    run(paths_file);
end
subj_dir = input_dir;   % paths.m provides input_dir/output_dir, not subj_dir
run(fullfile(mimm_root, 'MIMM_set_path.m'));

%% --- Sweep configuration ---
lambda_list   = [0, 0.005, 0.01, 0.015, 0.025, 0.05, 0.1, 0.2, 0.4];
lambda_default = 0.015;
TE = [6, 12, 18, 24, 30] * 1e-3;   % seconds (MUMC protocol)

% Lock the knee as the cohort-wide lambda_chi (shared by all subjects) when
% running on your chosen reference subject. The defensible workflow: derive
% once, apply fixed to everyone — avoids the circularity of re-tuning
% lambda_chi on each subject's own QSM.
%   To lock without editing this file:  lock_cohort = true; run('sweep_lambda_chi.m')
if ~exist('lock_cohort', 'var'); lock_cohort = false; end

%% --- Load dictionary (MUMC if present, else original) ---
mumc_dict = fullfile(mimm_root, 'Dictionary', 'MIMM_dictionary_stochastic_MUMC.mat');
orig_dict = fullfile(mimm_root, 'Dictionary', 'MIMM_dictionary_stochastic.mat');
if exist(mumc_dict, 'file')
    stoch = load(mumc_dict, 'dictionary');
    disp('Using MUMC dictionary.');
else
    stoch = load(orig_dict, 'dictionary');
    warning('MUMC dictionary not found — using original with TE interpolation.');
end
dictionary = stoch.dictionary;

%% --- Load data ---
mag_data   = double(niftiread(fullfile(subj_dir, 'magnitude.nii.gz')));
QSM        = double(niftiread(fullfile(subj_dir, 'qsm', 'QSM.nii.gz')));
Brain_Mask = logical(niftiread(fullfile(subj_dir, 'qsm', 'brain_mask.nii.gz')));

%% --- Dictionary preparation (mirrors MIMM.m lines 36-59) ---
dictionary.magnitude = abs(dictionary.signal);
TE_original = dictionary.TE * 1e3;   % ms
TE_target   = TE * 1e3;              % ms
degree = 5;
mag_interp = interpolate_dictionary(dictionary.magnitude, TE_original, TE_target, degree);
dict_norm  = mag_interp ./ vecnorm(mag_interp, 2, 2);   % [Ndict x Necho]

chi_iron = dictionary.chi_iron * 1e6;   % ppm
chi_iso  = dictionary.chi_iso  * 1e6;   % ppm
dict_MVF = dictionary.FVF .* (1 - dictionary.g_ratio.^2);
chi_neg  = chi_iso .* dict_MVF;          % diamagnetic (negative)
chi_pos  = dictionary.IVF * chi_iron;    % paramagnetic (positive)
chi_total_dict = (chi_neg + chi_pos).';  % [Ndict x 1]

%% --- Data preparation (mirrors MIMM.m lines 65-85) ---
magnitude = abs(mag_data) .* Brain_Mask;
magnitude = magnitude ./ vecnorm(magnitude, 2, 4);
magnitude(isnan(magnitude)) = 0;
magnitude(isinf(magnitude)) = 0;

lims = bounding_box(Brain_Mask);
BM   = Brain_Mask(lims{:});
QSMc = QSM(lims{:});
magc = magnitude(lims{:}, :) .* BM;
s0   = size(magc);

K = numel(lambda_list);

% Accumulate across all in-brain voxels, per lambda.
all_qsm = {};   % collected once (lambda-independent)
chi_est_per_lambda = cell(1, K);   % chi_total estimate
magfit_per_lambda  = cell(1, K);   % achieved magnitude-decay correlation
for k = 1:K; chi_est_per_lambda{k} = {}; magfit_per_lambda{k} = {}; end

%% --- Matching pursuit: compute error terms once per slice, loop lambda ---
fprintf('Sweeping %d lambda values over %d slices...\n', K, s0(3));
for slice = 1:s0(3)
    bm_slice = reshape(BM(:,:,slice), prod(s0(1:2)), 1);
    vox = find(bm_slice);            % in-brain voxels only
    if isempty(vox); continue; end

    mag_slice = reshape(magc(:,:,slice,:), prod(s0(1:2)), s0(4)).';
    qsm_slice = reshape(QSMc(:,:,slice), prod(s0(1:2)), 1).';
    mag_slice = mag_slice(:, vox);   % [Necho x Nvox_mask]
    qsm_slice = qsm_slice(vox);      % [1 x Nvox_mask]

    % Lambda-independent error terms (the expensive part).
    mag_corr  = abs(dict_norm * mag_slice);           % [Ndict x Nvox]
    base_err  = 1 - mag_corr;                          % [Ndict x Nvox]
    chi_err   = abs(chi_total_dict - qsm_slice);       % [Ndict x Nvox]

    all_qsm{end+1} = qsm_slice(:);                     %#ok<SAGROW>
    nv = numel(vox);

    % Cheap argmin per lambda; record chi estimate and achieved magnitude fit.
    for k = 1:K
        [~, ind] = min(lambda_list(k) * chi_err + base_err, [], 1);
        chi_est_per_lambda{k}{end+1} = chi_total_dict(ind);
        lin = sub2ind(size(mag_corr), ind, 1:nv);
        magfit_per_lambda{k}{end+1}  = mag_corr(lin).';
    end

    if mod(slice, 20) == 0
        fprintf('  slice %d/%d\n', slice, s0(3));
    end
end

qsm_vec = cell2mat(all_qsm(:));

%% --- Per-lambda statistics ---
bias  = zeros(K,1); slope = zeros(K,1); rval = zeros(K,1);
mae   = zeros(K,1); magfit = zeros(K,1);
for k = 1:K
    est = cell2mat(chi_est_per_lambda{k}(:));
    d   = est - qsm_vec;
    bias(k)  = mean(d);
    mae(k)   = mean(abs(d));
    p        = polyfit(qsm_vec, est, 1);
    slope(k) = p(1);
    R        = corrcoef(qsm_vec, est);
    rval(k)  = R(1,2);
    magfit(k)= mean(cell2mat(magfit_per_lambda{k}(:)));   % mean magnitude-decay fit
    fprintf('lambda=%.3f : bias=%+.4f  slope=%.3f  r=%.3f  MAE=%.4f  magfit=%.4f\n', ...
        lambda_list(k), bias(k), slope(k), rval(k), mae(k), magfit(k));
end

%% --- Detect L-curve knee (utopia-point method) ---
% Normalise both axes to [0,1]; the ideal corner is (low chi MAE, high mag fit).
% The knee is the sampled lambda closest to that corner = best trade-off.
mae_n    = (mae    - min(mae))    ./ max(max(mae)    - min(mae),    eps);
magfit_n = (magfit - min(magfit)) ./ max(max(magfit) - min(magfit), eps);
dist2corner = sqrt(mae_n.^2 + (1 - magfit_n).^2);
[~, knee_idx] = min(dist2corner);
knee_lambda   = lambda_list(knee_idx);
fprintf('\nL-curve knee: lambda_chi = %.4g  (slope=%.3f, bias=%+.4f, magfit=%.4f)\n', ...
    knee_lambda, slope(knee_idx), bias(knee_idx), magfit(knee_idx));

% Write the per-subject recommendation as a data file. run_MIMM_MUMC.m reads
% this if no cohort lock exists; it does NOT overwrite the source default.
rec_file = fullfile(mimm_dir, 'lambda_chi_recommended.txt');
fid = fopen(rec_file, 'w');
fprintf(fid, '%.6f\n', knee_lambda);
fclose(fid);
fprintf('Saved: %s\n', rec_file);

% If locking, also write the cohort-wide value next to paths.m, where every
% subject's run_MIMM_MUMC.m can find it. This takes precedence over per-subject.
if lock_cohort
    cohort_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'lambda_chi_cohort.txt');
    fid = fopen(cohort_file, 'w');
    fprintf(fid, '%.6f\n', knee_lambda);
    fclose(fid);
    fprintf('LOCKED cohort lambda_chi = %.4g -> %s\n', knee_lambda, cohort_file);
    fprintf('  All subjects will now use this value. Delete the file to unlock.\n');
end

%% --- Save CSV ---
out_csv = fullfile(mimm_dir, 'lambda_chi_sweep.csv');
fid = fopen(out_csv, 'w');
fprintf(fid, 'lambda_chi,bias,slope,r,mean_abs_err,magnitude_fit\n');
for k = 1:K
    fprintf(fid, '%.4f,%.6f,%.6f,%.6f,%.6f,%.6f\n', ...
        lambda_list(k), bias(k), slope(k), rval(k), mae(k), magfit(k));
end
fclose(fid);
fprintf('Saved: %s\n', out_csv);

%% --- Plot: slope / bias / magnitude fit vs lambda + L-curve ---
fig = figure('Visible','off','Position',[0 0 1100 850],'Color','w');

subplot(2,2,1);
semilogx(lambda_list, slope, '-o', 'LineWidth', 1.5); hold on;
yline(1, '--', 'ideal (slope=1)');
xline(lambda_default, ':', 'default');
xlabel('\lambda_{chi}'); ylabel('slope (\chi_{total} vs QSM)');
title('Slope \rightarrow 1 = QSM range reproduced'); grid on;

subplot(2,2,2);
semilogx(lambda_list, bias, '-o', 'LineWidth', 1.5); hold on;
yline(0, '--', 'ideal (bias=0)');
xline(lambda_default, ':', 'default');
xlabel('\lambda_{chi}'); ylabel('bias (\chi_{total} - QSM) [ppm]');
title('Bias \rightarrow 0 = no offset'); grid on;

subplot(2,2,3);
semilogx(lambda_list, magfit, '-o', 'LineWidth', 1.5); hold on;
xline(lambda_default, ':', 'default');
xlabel('\lambda_{chi}'); ylabel('mean magnitude-decay fit');
title('Magnitude fit (drops as \lambda forces QSM)'); grid on;

% L-curve: magnitude fit (quality we sacrifice) vs chi MAE (quality we gain).
% The knee = best trade-off = data-driven optimal lambda_chi.
subplot(2,2,4);
plot(mae, magfit, '-o', 'LineWidth', 1.5); hold on;
plot(mae(knee_idx), magfit(knee_idx), 'p', 'MarkerSize', 16, ...
     'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'k');   % knee = red star
for k = 1:K
    if lambda_list(k) == lambda_default
        txt = sprintf(' %.3g (default)', lambda_list(k));
    else
        txt = sprintf(' %.3g', lambda_list(k));
    end
    text(mae(k), magfit(k), txt, 'FontSize', 8);
end
xlabel('\chi MAE (|\chi_{total} - QSM|) [ppm]'); ylabel('mean magnitude-decay fit');
title(sprintf('L-curve: knee at \\lambda_{chi} = %.3g', knee_lambda));
grid on; set(gca,'XDir','reverse');

sgtitle('MIMM \lambda_{chi} sweep: \chi_{total} vs measured QSM');
out_png = fullfile(mimm_dir, 'lambda_chi_sweep.png');
saveas(fig, out_png);
close(fig);
fprintf('Saved: %s\n', out_png);

fprintf(['\nReading the result:\n' ...
    '  - If slope rises toward 1 and bias toward 0 as lambda_chi increases,\n' ...
    '    the default 0.015 is too weak and a larger value fits QSM better.\n' ...
    '  - Magnitude fit (panel 3) drops as lambda forces QSM agreement; this is\n' ...
    '    the cost (MVF accuracy depends on the decay-curve match).\n' ...
    '  - L-curve (panel 4): the knee is the data-driven optimal lambda_chi —\n' ...
    '    maximal QSM agreement before magnitude fit collapses.\n']);
