% Extract subject-level ROI means from MIMM, MWF, and chi-separation maps
%
% For each subject: loads output NIfTI maps, applies ROI masks, computes
% the mean of each parameter within each ROI. Produces [N_subjects x 1]
% vectors — the correct input for Bland-Altman and clinical correlation.
%
% ROI masks expected per subject (adapt to your atlas/segmentation pipeline):
%   roi/NAWM_mask.nii.gz       — normal-appearing white matter
%   roi/lesion_mask.nii.gz     — white matter lesions (MS subjects only)
%
% All maps must already be registered to the same space (ME-GRE 1 mm iso).

%% --- Study paths (adapt) ---

subj_list  = {             % list all subject directories
    '/data/MUMC/sub-01',
    '/data/MUMC/sub-02',
    % ...
};

output_file = 'roi_means.mat';   % saved in current directory

%% --- Configuration ---

mimm_variant = 'basic';   % 'basic', 'DTI', or 'Atlas'

roi_names = {'NAWM', 'lesion'};
roi_files = {'roi/NAWM_mask.nii.gz', 'roi/lesion_mask.nii.gz'};

mimm_fields  = {'MVF', 'FVF', 'g_ratio', 'R2s', 'chi_iron_est', 'chi_myelin'};
other_fields = {'MWF', 'chi_neg'};   % from T2-GRASE and chi-separation

%% --- Initialise output struct ---

N = numel(subj_list);
for ri = 1:numel(roi_names)
    roi = roi_names{ri};
    for fi = 1:numel(mimm_fields)
        roi_means.(roi).(mimm_fields{fi}) = nan(N, 1);
    end
    for fi = 1:numel(other_fields)
        roi_means.(roi).(other_fields{fi}) = nan(N, 1);
    end
end

%% --- Loop over subjects ---

for si = 1:N
    sdir = subj_list{si};
    fprintf('Subject %d/%d: %s\n', si, N, sdir);

    % --- Load MIMM maps ---
    mimm_dir = fullfile(sdir, 'mimm');
    mimm_maps = struct();
    all_ok = true;
    for fi = 1:numel(mimm_fields)
        fn = fullfile(mimm_dir, sprintf('%s_%s.nii.gz', mimm_fields{fi}, mimm_variant));
        if ~exist(fn, 'file')
            warning('Missing MIMM map: %s', fn);
            all_ok = false; break;
        end
        mimm_maps.(mimm_fields{fi}) = double(niftiread(fn));
    end
    if ~all_ok; continue; end

    % --- Load MWF map (registered to ME-GRE space) ---
    mwf_file = fullfile(sdir, 't2grase', 'MWF_reg.nii.gz');
    if exist(mwf_file, 'file')
        MWF = double(niftiread(mwf_file));
    else
        warning('Missing MWF map: %s', mwf_file);
        MWF = [];
    end

    % --- Load chi_neg map (from chi-separation) ---
    chineg_file = fullfile(sdir, 'source_separation', 'chi_neg.nii.gz');
    if exist(chineg_file, 'file')
        chi_neg = double(niftiread(chineg_file));
    else
        warning('Missing chi_neg map: %s', chineg_file);
        chi_neg = [];
    end

    % --- Extract means per ROI ---
    for ri = 1:numel(roi_names)
        roi      = roi_names{ri};
        roi_file = fullfile(sdir, roi_files{ri});

        if ~exist(roi_file, 'file')
            warning('Missing ROI mask: %s', roi_file);
            continue;
        end

        mask = logical(niftiread(roi_file));
        if sum(mask(:)) == 0
            warning('Empty ROI mask: %s', roi_file);
            continue;
        end

        % MIMM parameters
        for fi = 1:numel(mimm_fields)
            vals = mimm_maps.(mimm_fields{fi})(mask);
            roi_means.(roi).(mimm_fields{fi})(si) = mean(vals(isfinite(vals)));
        end

        % MWF
        if ~isempty(MWF)
            vals = MWF(mask);
            roi_means.(roi).MWF(si) = mean(vals(isfinite(vals)));
        end

        % chi_neg (use absolute value — diamagnetic so negative by convention)
        if ~isempty(chi_neg)
            vals = abs(chi_neg(mask));
            roi_means.(roi).chi_neg(si) = mean(vals(isfinite(vals)));
        end
    end
end

%% --- Save ---

save(output_file, 'roi_means', 'subj_list', 'mimm_variant');
fprintf('ROI means saved to: %s\n', output_file);

%% --- Quick summary table ---

fprintf('\nSummary (mean ± SD across subjects)\n');
fprintf('%-12s %-12s %s\n', 'ROI', 'Parameter', 'Mean ± SD');
for ri = 1:numel(roi_names)
    roi = roi_names{ri};
    all_params = [mimm_fields, other_fields];
    for fi = 1:numel(all_params)
        p = all_params{fi};
        v = roi_means.(roi).(p);
        v = v(~isnan(v));
        if ~isempty(v)
            fprintf('%-12s %-12s %.4f ± %.4f  (n=%d)\n', roi, p, mean(v), std(v), numel(v));
        end
    end
end
