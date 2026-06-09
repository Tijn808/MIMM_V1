% QSM reconstruction for MUMC ME-GRE data using the Chi-separation toolbox pipeline.
% Uses ROMEO phase unwrapping + V-SHARP + iLSQR (Method 1 from Chisep_script.m).
%
% Run setup_toolbox_paths.m once before this script.
%
% Inputs (in subj_dir):
%   magnitude.nii.gz   [224 x 224 x 154 x 5]   (prepared by prepare_mgre.m)
%   phase.nii.gz       [224 x 224 x 154 x 5]   (radians, prepared by prepare_mgre.m)
%
% Outputs (written to subj_dir/qsm/):
%   QSM.nii.gz          susceptibility map (ppm)
%   R2star.nii.gz       R2* map (s^-1)
%   brain_mask.nii.gz   binary brain mask

%% --- Paths (loaded from MUMC_pipeline/paths.m) ---

if ~exist('mimm_root', 'var')   % skip if called from run_subject.m with vars pre-set
    paths_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'paths.m');
    if ~exist(paths_file, 'file')
        error('paths.m not found. Copy MUMC_pipeline/paths_template.m to paths.m and fill in your paths.');
    end
    run(paths_file);
end
subj_dir   = input_dir;
output_dir = qsm_dir;

%% --- Load toolboxes ---
% Prefer the toolbox's own setup script if present. Otherwise (e.g. a copy that
% does not ship setup_toolbox_paths.m, or one nested inside extra folders),
% configure it ourselves: add everything under chisep_dir, then make sure the
% ROMEO build whose binary actually runs on THIS machine is first on the path
% (the 22.04 build runs on glibc 2.34; the 24.04 one needs a newer glibc).
setup_file = fullfile(chisep_dir, 'setup_toolbox_paths.m');
if exist(setup_file, 'file')
    run(setup_file);
else
    fprintf('setup_toolbox_paths.m not found in %s — auto-configuring.\n', chisep_dir);
    addpath(genpath(chisep_dir));
    romeo_ms = dir(fullfile(chisep_dir, '**', 'mritools_*', 'matlab', 'ROMEO.m'));
    picked = '';
    for k = 1:numel(romeo_ms)
        rbin = fullfile(romeo_ms(k).folder, '..', 'bin', 'romeo');
        if isfile(rbin)
            [st, ~] = system(sprintf('"%s" --help', rbin));
            if st == 0
                addpath(romeo_ms(k).folder);   % prepend: this ROMEO.m wins
                picked = rbin;
                break;
            end
        end
    end
    if isempty(picked)
        error(['No working ROMEO binary found under %s. ' ...
               'Checked %d mritools build(s); none ran on this machine.'], ...
               chisep_dir, numel(romeo_ms));
    end
    fprintf('Using ROMEO binary: %s\n', picked);
end

%% --- Acquisition parameters ---

Data = struct();
Data.TE           = [6 12 18 24 30];    % ms
Data.B0dir        = [0 0 1];
Data.CF           = 127e6;              % Hz (3T: 42.577e6 * 3)
Data.B0_strength  = 3;                  % T
Data.VoxelSize    = [1 1 1];            % mm

RunOptions = struct();
RunOptions.InputType    = 'nifti';
RunOptions.Mask         = false;
RunOptions.Mask_method  = 'MEDI';
RunOptions.R2sfit       = 'ARLO';
RunOptions.Unwrap       = 'ROMEO + weighted echo averaging';
RunOptions.BFR          = 'V-SHARP';
RunOptions.Tukey        = 0;            % Philips data: no Tukey windowing
RunOptions.PhaseInverse = 0;
RunOptions.denoising    = false;
RunOptions.EvenSizePadding = [0 0 0];
RunOptions.OutputPath   = output_dir;

Data.RunOptions = RunOptions;
Data.output_root = output_dir;

%% --- Load NIfTI data ---

disp('Loading magnitude and phase...')
nii_mag = niftiread(fullfile(subj_dir, 'magnitude.nii.gz'));
nii_phs = niftiread(fullfile(subj_dir, 'phase.nii.gz'));

Data.MGRE_Mag = double(nii_mag);   % [x y z echoes], already positive
Data.MGRE_Phs = double(nii_phs);   % radians [-pi, pi] from prepare_mgre.m

Data.MatrixSize = size(Data.MGRE_Mag);
Data.Necho      = size(Data.MGRE_Mag, 4);

% Force even dimensions (required by chi-sep internals)
input_fields = {'MGRE_Mag', 'MGRE_Phs'};
for i = 1:length(input_fields)
    [Data.(input_fields{i}), x_odd, y_odd, z_odd] = even_pad(Data.(input_fields{i}));
end
RunOptions.EvenSizePadding = [x_odd, y_odd, z_odd];
Data.MatrixSize = size(Data.MGRE_Mag);

% TE shape: must be column vector
Data.TE = Data.TE(:);

%% --- Tukey windowing (disabled for Philips, Tukey=0 = no windowing) ---

% Tukey=0 means no windowing — skip tukey_windowing (requires Signal Processing Toolbox)
Data.MGRE_Mag_Tukey = Data.MGRE_Mag;
Data.MGRE_Phs_Tukey = Data.MGRE_Phs * (-1)^RunOptions.PhaseInverse;

%% --- Brain mask ---

disp('Running brain extraction (MEDI BET)...')
Data.Mask = double(BET(Data.MGRE_Mag_Tukey(:,:,:,1), Data.MatrixSize(1:3), Data.VoxelSize));

%% --- R2* fitting (ARLO) ---

disp('Fitting R2* (ARLO)...')
Data.dTE = Data.TE(2) - Data.TE(1);   % ms
Data.R2s = r2star_arlo(Data.MGRE_Mag_Tukey, Data.TE, Data.Mask);
Data.R2s(Data.R2s < 0) = 0;
Data.map = Data.R2s;

%% --- Phase unwrapping (ROMEO) ---

disp('Unwrapping phase with ROMEO...')
parameters.TE                    = Data.TE';
parameters.mag                   = Data.MGRE_Mag_Tukey;
parameters.mask                  = double(Data.Mask);
parameters.calculate_B0          = false;
parameters.phase_offset_correction = 'on';
parameters.voxel_size            = Data.VoxelSize;
parameters.additional_flags      = '-q';
parameters.output_dir            = fullfile(output_dir, 'romeo_tmp');
mkdir(parameters.output_dir);

[unwrapped_phase, ~] = ROMEO(double(Data.MGRE_Phs_Tukey), parameters);
unwrapped_phase(isnan(unwrapped_phase)) = 0;

% Weighted echo averaging
TE_s     = Data.TE / 1000;
t2s_roi  = 0.04;
W        = TE_s .* exp(-TE_s / t2s_roi);
W        = W / sum(W);
weighted = zeros(Data.MatrixSize(1:3));
TE_eff   = 0;
for e = 1:Data.Necho
    weighted = weighted + W(e) * unwrapped_phase(:,:,:,e);
    TE_eff   = TE_eff   + W(e) * TE_s(e);
end
Data.UnwrappedPhase = weighted / TE_eff * (TE_s(2) - TE_s(1)) .* Data.Mask;

%% --- Background field removal (V-SHARP) ---

disp('Removing background field (V-SHARP)...')
[Data.local_field, Data.mask_brain_new] = V_SHARP(Data.UnwrappedPhase, Data.Mask, ...
    'voxelsize', Data.VoxelSize, 'smvsize', 12);
Data.delta_TE      = (Data.TE(2) - Data.TE(1)) / 1000;   % seconds
Data.local_field_hz = double(Data.local_field) / (2 * pi * Data.delta_TE);

%% --- QSM via iLSQR (STI Suite) ---

disp('Running iLSQR dipole inversion...')
pad_size = [12 12 12];
% QSM_iLSQR expects local_field in radians (V_SHARP output), voxelsize as row vector
Data.QSM = QSM_iLSQR(Data.local_field, Data.mask_brain_new, ...
    'TE',       Data.delta_TE * 1e3, ...
    'B0',       Data.B0_strength, ...
    'H',        Data.B0dir, ...
    'padsize',  pad_size, ...
    'voxelsize', Data.VoxelSize);

%% --- Unpad if even-size padding was applied ---

if any(RunOptions.EvenSizePadding)
    fields_to_unpad = {'QSM', 'R2s', 'mask_brain_new'};
    for i = 1:length(fields_to_unpad)
        if isfield(Data, fields_to_unpad{i})
            Data.(fields_to_unpad{i}) = even_unpad(Data.(fields_to_unpad{i}), RunOptions.EvenSizePadding);
        end
    end
end

QSM        = Data.QSM;
R2star     = Data.R2s;
Brain_Mask = logical(Data.mask_brain_new);

%% --- Save outputs ---

% Use a clean dummy header (no inherited Philips scl_slope from magnitude)
ref_info = niftiinfo(fullfile(subj_dir, 'magnitude.nii.gz'));
dummy_file = fullfile(tempdir, 'chisep_ref3d.nii');
niftiwrite(single(QSM), dummy_file);
ref3d = niftiinfo(dummy_file);
delete(dummy_file);
ref3d.PixelDimensions = ref_info.PixelDimensions(1:3);
ref3d.SpaceUnits      = 'Millimeter';

ref_dbl      = ref3d; ref_dbl.Datatype = 'double'; ref_dbl.BitsPerPixel = 64;
ref_u8       = ref3d; ref_u8.Datatype  = 'uint8';  ref_u8.BitsPerPixel  = 8;

niftiwrite(double(QSM),          fullfile(output_dir, 'QSM.nii.gz'),        ref_dbl, 'Compressed', true);
niftiwrite(double(R2star),       fullfile(output_dir, 'R2star.nii.gz'),      ref_dbl, 'Compressed', true);
niftiwrite(uint8(Brain_Mask),    fullfile(output_dir, 'brain_mask.nii.gz'),  ref_u8,  'Compressed', true);

% Also save echo 1 magnitude for atlas registration
mag_e1 = double(niftiread(fullfile(subj_dir, 'magnitude.nii.gz')));
mag_e1 = mag_e1(:,:,:,1);
ref_e1 = ref3d; ref_e1.Datatype = 'double'; ref_e1.BitsPerPixel = 64;
niftiwrite(mag_e1, fullfile(output_dir, 'mag_e1.nii.gz'), ref_e1, 'Compressed', true);

% Clean up ROMEO temp files
rmdir(fullfile(output_dir, 'romeo_tmp'), 's');

mask_qs = logical(Brain_Mask);
qsm_vals = QSM(mask_qs);
r2s_vals = R2star(Brain_Mask);
fprintf('QSM  P5/P50/P95: [%.4f, %.4f, %.4f] ppm\n', prctile(qsm_vals,5), prctile(qsm_vals,50), prctile(qsm_vals,95));
fprintf('R2*  P5/P50/P95: [%.1f, %.1f, %.1f] s^-1\n', prctile(r2s_vals,5), prctile(r2s_vals,50), prctile(r2s_vals,95));
disp('QSM reconstruction complete.')
fprintf('Output written to: %s\n', output_dir)
