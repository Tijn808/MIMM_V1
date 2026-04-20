cd('/home/tijn-saes/Documents/Internship/MIMM')
MIMM_set_path

% Load results
results_stoch  = load('Example_Results/stochastic_MIMM_results_new.mat');
results_determ = load('Example_Results/deterministic_MIMM_results_new.mat');

datasets       = {results_stoch, results_determ};
dataset_labels = {'stochastic', 'deterministic'};

fields = {'MVF', 'g_ratio', 'FVF', 'R2s', 'chi_iron_est', 'chi_myelin', 'theta_est', 'error'};
clims  = {[0 0.4], [0.5 1], [0 0.8], [0 100], [0 0.3], [-0.06 0], [0 90], [0 0.25]};
modes       = {'MIMM_basic', 'MIMM_DTI', 'MIMM_Atlas'};
mode_labels = {'Basic', 'DTI', 'Atlas'};

for d = 1:length(datasets)
    results       = datasets{d};
    dataset_label = dataset_labels{d};

    output_dir = ['Example_Results/figures/' dataset_label];
    if ~exist(output_dir, 'dir')
        mkdir(output_dir)
    end

    slice = round(size(results.MIMM_basic.MVF, 3) / 2);

    for f = 1:length(fields)
        field = fields{f};
        fig = figure('Visible','off','Position',[0 0 1200 400]);

        for m = 1:length(modes)
            data = results.(modes{m}).(field)(:,:,slice);
            subplot(1,3,m)
            imagesc(data', clims{f})
            colormap(gca, 'gray')
            colorbar
            axis image off
            title(sprintf('%s - %s', mode_labels{m}, field), 'Interpreter','none')
        end

        saveas(fig, fullfile(output_dir, sprintf('%s_slice%d.png', field, slice)))
        close(fig)
        fprintf('Saved %s/%s\n', dataset_label, field)
    end
end

disp('All figures saved!')
