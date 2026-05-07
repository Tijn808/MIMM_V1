% Bland-Altman analysis: MIMM MVF vs MWF and vs source-separation myelin proxy
%
% IMPORTANT: inputs must be SUBJECT-LEVEL means (one value per subject per ROI).
% Do NOT pass voxel-level data — pooled voxels from multiple subjects violate
% the independence assumption, inflate N, and produce artificially narrow LoA.
%
% Inputs:
%   mvf_mimm    : [N_subjects x 1] mean MIMM MVF per subject within ROI
%   mwf         : [N_subjects x 1] mean T2-GRASE MWF per subject within same ROI
%   chi_neg     : [N_subjects x 1] mean |chi_neg| per subject within same ROI (ppm)
%   roi_label   : string label for the ROI (for figure titles)
%
% Outputs:
%   Bland-Altman plot saved as PNG
%   Struct with bias, LoA, and Pearson r

function ba = bland_altman(mvf_mimm, mwf, chi_neg, roi_label, output_dir)

    if nargin < 4; roi_label   = 'ROI'; end
    if nargin < 5; output_dir  = '.';   end

    % Enforce subject-level input (warn if suspiciously large N)
    if numel(mvf_mimm) > 200
        warning(['bland_altman received %d data points — expected N_subjects (typically <100).\n' ...
                 'Make sure you are passing subject-level ROI means, not voxel-level data.'], ...
                numel(mvf_mimm));
    end

    %% --- MIMM MVF vs MWF ---
    [ba.mvf_vs_mwf, fig1] = ba_pair(mvf_mimm, mwf, ...
        'MIMM MVF', 'T2-GRASE MWF', roi_label);
    saveas(fig1, fullfile(output_dir, ['BA_MVF_vs_MWF_' roi_label '.png']));
    close(fig1);

    %% --- MIMM MVF vs chi_neg (source separation) ---
    % chi_neg is in ppm; scale to comparable range via linear regression
    % before Bland-Altman (the two measures are not in the same units)
    [ba.mvf_vs_chisep, fig2] = ba_pair_scaled(mvf_mimm, chi_neg, ...
        'MIMM MVF', '|chi_{neg}| (ppm)', roi_label);
    saveas(fig2, fullfile(output_dir, ['BA_MVF_vs_ChiSep_' roi_label '.png']));
    close(fig2);

    disp(['Bland-Altman complete for ROI: ' roi_label]);
    disp_ba_stats(ba.mvf_vs_mwf,    'MVF vs MWF');
    disp_ba_stats(ba.mvf_vs_chisep, 'MVF vs |chi_neg|');
end


%% --- Helper: same-unit Bland-Altman ---
function [stats, fig] = ba_pair(x, y, xlabel_str, ylabel_str, roi_label)
    x = x(:); y = y(:);
    mean_xy = (x + y) / 2;
    diff_xy = x - y;

    bias = mean(diff_xy);
    sd   = std(diff_xy);
    loa_upper = bias + 1.96 * sd;
    loa_lower = bias - 1.96 * sd;

    % Pearson correlation
    r = corr(x, y, 'type', 'Pearson');

    stats.bias      = bias;
    stats.sd        = sd;
    stats.loa_upper = loa_upper;
    stats.loa_lower = loa_lower;
    stats.r         = r;
    stats.n         = numel(x);

    fig = figure('Visible', 'off');
    scatter(mean_xy, diff_xy, 6, [0.2 0.4 0.8], 'filled', 'MarkerFaceAlpha', 0.4);
    hold on;
    yline(bias,      'k-',  sprintf('Bias = %.3f', bias),      'LabelHorizontalAlignment', 'left');
    yline(loa_upper, 'r--', sprintf('+1.96 SD = %.3f', loa_upper), 'LabelHorizontalAlignment', 'left');
    yline(loa_lower, 'r--', sprintf('-1.96 SD = %.3f', loa_lower), 'LabelHorizontalAlignment', 'left');
    xlabel(['Mean of ' xlabel_str ' and ' ylabel_str]);
    ylabel([xlabel_str ' \minus ' ylabel_str]);
    title(sprintf('Bland-Altman: %s vs %s — %s\nr = %.3f, n = %d', ...
          xlabel_str, ylabel_str, roi_label, r, numel(x)));
    box on; grid on;
end


%% --- Helper: cross-unit Bland-Altman (scale y to x via OLS before differencing) ---
function [stats, fig] = ba_pair_scaled(x, chi_neg, xlabel_str, ylabel_str, roi_label)
    x       = x(:);
    chi_neg = abs(chi_neg(:));   % ensure positive

    % OLS: x = a + b * chi_neg  → scale chi_neg to MVF units
    b_coef  = [ones(size(chi_neg)), chi_neg] \ x;
    chi_scaled = b_coef(1) + b_coef(2) * chi_neg;

    [stats, fig] = ba_pair(x, chi_scaled, xlabel_str, [ylabel_str ' (scaled)'], roi_label);
    stats.scale_intercept = b_coef(1);
    stats.scale_slope     = b_coef(2);
end


%% --- Helper: print stats to console ---
function disp_ba_stats(stats, label)
    fprintf('%s: bias=%.4f, SD=%.4f, LoA=[%.4f, %.4f], r=%.3f, n=%d\n', ...
        label, stats.bias, stats.sd, stats.loa_lower, stats.loa_upper, stats.r, stats.n);
end
