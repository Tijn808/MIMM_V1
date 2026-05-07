% Spearman correlation: MIMM parameters vs clinical scores
%
% Inputs:
%   params   : struct with fields mvf, fvf, g_ratio, chi_iron — each [N_subjects x 1]
%              (mean value per subject within a chosen ROI, or whole-brain mean)
%   scores   : struct with fields edss, t25fw, nhpt, ssdmt — each [N_subjects x 1]
%              (NaN for missing values)
%   roi_label: string, used in figure titles and filenames
%   output_dir: directory to save PNG files
%
% Outputs:
%   results  : struct with rho and p-value for each param x score combination
%   Scatter plots and heatmap saved as PNG

function results = correlate_clinical(params, scores, roi_label, output_dir)

    if nargin < 3; roi_label   = 'ROI'; end
    if nargin < 4; output_dir  = '.';   end

    param_names = {'mvf',     'fvf',    'g_ratio', 'chi_iron'};
    score_names = {'edss',    't25fw',  'nhpt',    'ssdmt'};
    score_labels = {'EDSS',  'T25FW (s)', '9HPT (s)', 'sSDMT (score)'};
    param_labels = {'MVF',   'FVF',    'g-ratio',  '\chi_{iron} (ppm)'};

    n_params = numel(param_names);
    n_scores = numel(score_names);

    rho_mat = nan(n_params, n_scores);
    p_mat   = nan(n_params, n_scores);

    %% --- Pairwise Spearman, individual scatter plots ---
    for pi = 1:n_params
        x = params.(param_names{pi});
        if isempty(x) || all(isnan(x)); continue; end

        for si = 1:n_scores
            y = scores.(score_names{si});
            if isempty(y) || all(isnan(y)); continue; end

            % Remove subjects with NaN in either variable
            valid = ~isnan(x) & ~isnan(y);
            if sum(valid) < 5
                warning('Too few valid pairs for %s vs %s (%d subjects)', ...
                    param_names{pi}, score_names{si}, sum(valid));
                continue;
            end

            [rho, p] = corr(x(valid), y(valid), 'type', 'Spearman');
            rho_mat(pi, si) = rho;
            p_mat(pi, si)   = p;

            % Scatter
            fig = figure('Visible', 'off');
            scatter(x(valid), y(valid), 40, [0.2 0.5 0.8], 'filled', 'MarkerFaceAlpha', 0.7);
            xlabel([param_labels{pi} ' (' roi_label ')']);
            ylabel(score_labels{si});
            title(sprintf('%s vs %s — %s\n\\rho = %.3f, p = %.4f, n = %d', ...
                param_labels{pi}, score_labels{si}, roi_label, rho, p, sum(valid)));
            add_regression_line(x(valid), y(valid));
            box on; grid on;
            fname = sprintf('Corr_%s_vs_%s_%s.png', param_names{pi}, score_names{si}, roi_label);
            saveas(fig, fullfile(output_dir, fname));
            close(fig);
        end
    end

    %% --- Heatmap of all rho values ---
    fig_hm = figure('Visible', 'off', 'Position', [100 100 600 400]);
    imagesc(rho_mat, [-1 1]);
    colormap(redblue_colormap());
    colorbar;
    set(gca, 'XTick', 1:n_scores, 'XTickLabel', score_labels, ...
             'YTick', 1:n_params, 'YTickLabel', param_labels);
    xtickangle(30);
    title(sprintf('Spearman \\rho — %s', roi_label));

    % Annotate cells with rho (bold if significant)
    for pi = 1:n_params
        for si = 1:n_scores
            if ~isnan(rho_mat(pi, si))
                if p_mat(pi, si) < 0.05
                    txt = sprintf('%.2f*', rho_mat(pi, si));
                    fw  = 'bold';
                else
                    txt = sprintf('%.2f',  rho_mat(pi, si));
                    fw  = 'normal';
                end
                text(si, pi, txt, 'HorizontalAlignment', 'center', ...
                     'FontWeight', fw, 'FontSize', 9, 'Color', 'k');
            end
        end
    end

    saveas(fig_hm, fullfile(output_dir, ['Heatmap_Spearman_' roi_label '.png']));
    close(fig_hm);

    %% --- Print summary table ---
    fprintf('\nSpearman correlations — %s\n', roi_label);
    fprintf('%-12s', '');
    for si = 1:n_scores; fprintf('%-14s', score_labels{si}); end
    fprintf('\n');
    for pi = 1:n_params
        fprintf('%-12s', param_labels{pi});
        for si = 1:n_scores
            if isnan(rho_mat(pi, si))
                fprintf('%-14s', 'n/a');
            else
                sig = '';
                if p_mat(pi, si) < 0.05;  sig = '*';  end
                if p_mat(pi, si) < 0.01;  sig = '**'; end
                if p_mat(pi, si) < 0.001; sig = '***'; end
                fprintf('%-14s', sprintf('%.3f%s (p=%.3f)', rho_mat(pi,si), sig, p_mat(pi,si)));
            end
        end
        fprintf('\n');
    end

    results.rho       = rho_mat;
    results.p         = p_mat;
    results.param_names = param_names;
    results.score_names = score_names;
end


%% --- Helpers ---

function add_regression_line(x, y)
    p = polyfit(x, y, 1);
    xl = xlim;
    xfit = linspace(xl(1), xl(2), 100);
    hold on;
    plot(xfit, polyval(p, xfit), 'k-', 'LineWidth', 1.2);
end

function cmap = redblue_colormap()
    % Simple diverging red-white-blue colormap (64 levels)
    n = 32;
    r = [linspace(0.8, 1, n), linspace(1, 0.2, n)]';
    g = [linspace(0.2, 1, n), linspace(1, 0.2, n)]';
    b = [linspace(0.2, 1, n), linspace(1, 0.8, n)]';
    cmap = [r, g, b];
end
