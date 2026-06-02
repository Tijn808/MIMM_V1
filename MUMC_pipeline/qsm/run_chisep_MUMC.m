% Chi-separation for MUMC ME-GRE data.
% Produces chi_pos (paramagnetic, iron) and chi_neg (diamagnetic, myelin).
% Uses Chi-separation (MEDI) method — no pre-computed QSM needed.
%
% Pipeline: MEDI BET → ARLO R2* → ROMEO unwrapping → V-SHARP → chi-sep
%
% Prerequisites: run setup_toolbox_paths.m (adds MEDI, STI Suite, ROMEO, chi-sep)
%
% Inputs (in subj_dir/qsm/):
%   magnitude.nii.gz   [224 x 224 x 154 x 5]
%   phase.nii.gz       [224 x 224 x 154 x 5]  radians
%
% Outputs (written to subj_dir/chisep/):
%   chi_pos.nii.gz     paramagnetic susceptibility (ppm)
%   chi_neg.nii.gz     diamagnetic susceptibility  (ppm)
%   chi_tot.nii.gz     total susceptibility = QSM  (ppm)
%   R2s.nii.gz         R2* map (s^-1)
%   brain_mask.nii.gz  binary brain mask

%% --- Paths (loaded from MUMC_pipeline/paths.m) ---

if ~exist('mimm_root', 'var')
    paths_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'paths.m');
    if ~exist(paths_file, 'file')
        error('paths.m not found. Copy MUMC_pipeline/paths_template.m to paths.m and fill in your paths.');
    end
    run(paths_file);
end
subj_dir   = input_dir;
output_dir = chisep_out;

%% --- Toolboxes ---

run(fullfile(chisep_dir, 'setup_toolbox_paths.m'));

%% --- Acquisition parameters ---

TE          = [6 12 18 24 30];   % ms (column vector enforced below)
B0dir       = [0 0 1];
CF          = 127e6;             % Hz (3T)
B0_strength = 3;                 % T
VoxelSize   = [1 1 1];          % mm

%% --- Load data ---

disp('Loading magnitude and phase...')
MGRE_Mag = double(niftiread(fullfile(subj_dir, 'magnitude.nii.gz')));
MGRE_Phs = double(niftiread(fullfile(subj_dir, 'phase.nii.gz')));   % radians

MatrixSize = size(MGRE_Mag);
Necho      = size(MGRE_Mag, 4);
TE         = TE(:);   % column vector

% No Tukey windowing for Philips data
MGRE_Mag_Tukey = MGRE_Mag;
MGRE_Phs_Tukey = MGRE_Phs;

%% --- Force even dimensions ---

[MGRE_Mag_Tukey, x_odd, y_odd, z_odd] = even_pad(MGRE_Mag_Tukey);
[MGRE_Phs_Tukey, ~, ~, ~]             = even_pad(MGRE_Phs_Tukey);
EvenSizePadding = [x_odd, y_odd, z_odd];
MatrixSize = size(MGRE_Mag_Tukey);

%% --- Brain mask (MEDI BET) ---

disp('Running brain extraction (MEDI BET)...')
Mask = double(BET(MGRE_Mag_Tukey(:,:,:,1), MatrixSize(1:3), VoxelSize));

%% --- R2* fitting (ARLO) ---

disp('Fitting R2* (ARLO)...')
dTE = TE(2) - TE(1);   % ms
R2s = r2star_arlo(MGRE_Mag_Tukey, TE, Mask);
R2s(R2s < 0) = 0;

%% --- Phase unwrapping (ROMEO) ---

disp('Unwrapping phase (ROMEO)...')
romeo_tmp = fullfile(output_dir, 'romeo_tmp');
mkdir(romeo_tmp);

parameters.TE                      = TE';
parameters.mag                     = MGRE_Mag_Tukey;
parameters.mask                    = double(Mask);
parameters.calculate_B0            = false;
parameters.phase_offset_correction = 'on';
parameters.voxel_size              = VoxelSize;
parameters.additional_flags        = '-q';
parameters.output_dir              = romeo_tmp;

[unwrapped_phase, ~] = ROMEO(double(MGRE_Phs_Tukey), parameters);
unwrapped_phase(isnan(unwrapped_phase)) = 0;

% Weighted echo averaging
TE_s    = TE / 1000;
t2s_roi = 0.04;
W       = TE_s .* exp(-TE_s / t2s_roi);
W       = W / sum(W);
weighted = zeros(MatrixSize(1:3));
TE_eff   = 0;
for e = 1:Necho
    weighted = weighted + W(e) * unwrapped_phase(:,:,:,e);
    TE_eff   = TE_eff   + W(e) * TE_s(e);
end
UnwrappedPhase = weighted / TE_eff * (TE_s(2) - TE_s(1)) .* Mask;

rmdir(romeo_tmp, 's');

%% --- Background field removal (V-SHARP) ---

disp('Background field removal (V-SHARP)...')
[local_field, mask_brain_new] = V_SHARP(UnwrappedPhase, Mask, ...
    'voxelsize', VoxelSize, 'smvsize', 12);
delta_TE       = (TE(2) - TE(1)) / 1000;   % seconds
local_field_hz = double(local_field) / (2 * pi * delta_TE);

%% --- CSF mask (needed by chi-sep MEDI) ---

mask_CSF = extract_CSF(R2s, mask_brain_new, VoxelSize);

%% --- Chi-separation (MEDI) ---

disp('Running chi-separation (MEDI)...')

mag_combined = sqrt(sum(MGRE_Mag_Tukey.^2, 4)) .* mask_brain_new;
lf_masked    = local_field_hz .* mask_brain_new;
r2prime      = R2s .* mask_brain_new;   % use R2* as proxy for R2'

[~, N_std] = Preprocessing4Phase(MGRE_Mag_Tukey, MGRE_Phs_Tukey);

params.b0_dir    = B0dir;
params.CF        = CF;
params.voxel_size = VoxelSize;
params.TE        = TE;
params.lambda    = 1;
params.lambda_CSF = 1;
params.Dr        = 137;

[x_para, x_dia, x_tot] = chi_sep_MEDI(mag_combined, lf_masked, r2prime, ...
    N_std, mask_brain_new, params, []);

%% --- Unpad if needed ---

if any(EvenSizePadding)
    x_para        = even_unpad(x_para,        EvenSizePadding);
    x_dia         = even_unpad(x_dia,         EvenSizePadding);
    x_tot         = even_unpad(x_tot,         EvenSizePadding);
    R2s           = even_unpad(R2s,           EvenSizePadding);
    mask_brain_new = even_unpad(mask_brain_new, EvenSizePadding);
end

%% --- Save outputs ---

disp('Saving outputs...')

% Use a clean 3D header (avoids inheriting 4D/scl_slope from magnitude)
dummy_file = fullfile(tempdir, 'chisep_ref.nii');
niftiwrite(single(x_para), dummy_file);
ref = niftiinfo(dummy_file);
delete(dummy_file);
ref.PixelDimensions = VoxelSize;
ref.SpaceUnits      = 'Millimeter';

ref_dbl      = ref; ref_dbl.Datatype = 'double'; ref_dbl.BitsPerPixel = 64;
ref_u8       = ref; ref_u8.Datatype  = 'uint8';  ref_u8.BitsPerPixel  = 8;

niftiwrite(double(x_para),         fullfile(output_dir, 'chi_pos.nii.gz'),    ref_dbl, 'Compressed', true);
niftiwrite(double(x_dia),          fullfile(output_dir, 'chi_neg.nii.gz'),    ref_dbl, 'Compressed', true);
niftiwrite(double(x_tot),          fullfile(output_dir, 'chi_tot.nii.gz'),    ref_dbl, 'Compressed', true);
niftiwrite(double(R2s),            fullfile(output_dir, 'R2s.nii.gz'),        ref_dbl, 'Compressed', true);
niftiwrite(uint8(mask_brain_new),  fullfile(output_dir, 'brain_mask.nii.gz'), ref_u8,  'Compressed', true);

mask_qs = logical(mask_brain_new);
fprintf('chi_pos P5/P50/P95: [%.4f  %.4f  %.4f] ppm\n', prctile(x_para(mask_qs),5), prctile(x_para(mask_qs),50), prctile(x_para(mask_qs),95));
fprintf('chi_neg P5/P50/P95: [%.4f  %.4f  %.4f] ppm\n', prctile(x_dia(mask_qs),5),  prctile(x_dia(mask_qs),50),  prctile(x_dia(mask_qs),95));
fprintf('R2*     P5/P50/P95: [%.1f  %.1f  %.1f] s^-1\n', prctile(R2s(mask_qs),5),    prctile(R2s(mask_qs),50),    prctile(R2s(mask_qs),95));
disp(['Chi-separation complete. Outputs in: ' output_dir])
