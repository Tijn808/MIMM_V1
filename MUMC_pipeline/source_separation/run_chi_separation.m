% Chi-separation: decompose QSM into paramagnetic (iron) and diamagnetic (myelin)
% susceptibility components.
%
% Method: Shin et al. 2021, NeuroImage 240:118371
%   R2'    = R2* - R2          (reversible transverse relaxation rate)
%   chi_pos = R2' / alpha_pos  (paramagnetic iron susceptibility)
%   chi_neg = QSM - chi_pos    (diamagnetic myelin susceptibility)
%
% Inputs (all registered to ME-GRE 1 mm iso space):
%   QSM map       : from qsm/run_QSM.m          (ppm)
%   R2* map       : from SEPIA total field step  (s^-1)
%   T2 map        : from T2-GRASE pipeline       (s)  — registered + resampled
%   Brain mask    : from qsm/brain_mask.nii.gz   (binary)
%
% Prerequisites:
%   T2-GRASE T2 map must be registered to ME-GRE space first.
%   Run in terminal (FSL flirt):
%     flirt -in t2grase/T2map.nii.gz \
%           -ref qsm/QSM.nii.gz \
%           -out source_separation/T2map_reg.nii.gz \
%           -omat source_separation/T2map_reg.mat \
%           -dof 6 -interp trilinear
%
% Output:
%   chi_pos.nii.gz  — paramagnetic susceptibility (iron, ppm, >= 0)
%   chi_neg.nii.gz  — diamagnetic susceptibility  (myelin, ppm, <= 0)
%   R2prime.nii.gz  — reversible relaxation rate  (s^-1, >= 0)

%% --- Parameters ---

% Iron relaxivity at 3T (s^-1 / ppm)
% From Shin et al. 2021: alpha_pos = 124 s^-1/ppm at 3T
alpha_pos = 124;

%% --- Paths (adapt per subject) ---

subj_dir   = '';    % e.g. '/data/MUMC/sub-01/'
output_dir = fullfile(subj_dir, 'source_separation');
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

%% --- Load inputs ---

qsm_file   = fullfile(subj_dir, 'qsm',              'QSM.nii.gz');
r2s_file   = fullfile(subj_dir, 'qsm',              'R2star.nii.gz');
t2_file    = fullfile(output_dir,                   'T2map_reg.nii.gz');
mask_file  = fullfile(subj_dir, 'qsm',              'brain_mask.nii.gz');

for f = {qsm_file, r2s_file, t2_file, mask_file}
    if ~exist(f{1}, 'file')
        error('Missing input: %s', f{1});
    end
end

QSM        = double(niftiread(qsm_file));
R2star     = double(niftiread(r2s_file));
T2         = double(niftiread(t2_file));
Brain_Mask = logical(niftiread(mask_file));

ref_info   = niftiinfo(qsm_file);

%% --- Compute R2 and R2' ---

% R2 = 1/T2 (avoid division by zero outside brain)
R2 = zeros(size(T2));
valid = T2 > 0 & Brain_Mask;
R2(valid) = 1 ./ T2(valid);

% R2' = R2* - R2, clamped to >= 0 (cannot be negative by definition)
R2prime = max(R2star - R2, 0) .* Brain_Mask;

%% --- Chi-separation ---

% Paramagnetic (iron): proportional to R2'
chi_pos = (R2prime / alpha_pos) .* Brain_Mask;

% Diamagnetic (myelin): remainder after subtracting iron contribution
chi_neg = (QSM - chi_pos) .* Brain_Mask;

% Convention: chi_neg should be <= 0 in myelin-containing tissue.
% Positive residuals in chi_neg indicate the iron model underestimates iron.
% Do NOT hard-clip chi_neg — preserve the full distribution for analysis.

%% --- Quality checks ---

fprintf('R2prime : mean=%.2f, max=%.2f s^-1 (brain)\n', ...
    mean(R2prime(Brain_Mask)), max(R2prime(Brain_Mask)));
fprintf('chi_pos : mean=%.4f, max=%.4f ppm (brain)\n', ...
    mean(chi_pos(Brain_Mask)), max(chi_pos(Brain_Mask)));
fprintf('chi_neg : mean=%.4f, min=%.4f ppm (brain)\n', ...
    mean(chi_neg(Brain_Mask)), min(chi_neg(Brain_Mask)));

pct_pos_chineg = 100 * sum(chi_neg(Brain_Mask) > 0) / sum(Brain_Mask(:));
fprintf('chi_neg > 0 in %.1f%% of brain voxels (expected ~0%% in pure WM)\n', pct_pos_chineg);

%% --- Save outputs ---

write_nii(R2prime, 'R2prime',  ref_info, output_dir);
write_nii(chi_pos, 'chi_pos',  ref_info, output_dir);
write_nii(chi_neg, 'chi_neg',  ref_info, output_dir);

disp('Chi-separation complete.');
fprintf('Outputs written to: %s\n', output_dir);


%% --- Helper ---
function write_nii(data, name, ref_info, out_dir)
    info = ref_info;
    info.Datatype      = 'double';
    info.BitsPerPixel  = 64;
    info.ImageSize     = size(data);
    niftiwrite(data, fullfile(out_dir, [name '.nii.gz']), info, 'Compressed', true);
    fprintf('  Written: %s.nii.gz\n', name);
end
