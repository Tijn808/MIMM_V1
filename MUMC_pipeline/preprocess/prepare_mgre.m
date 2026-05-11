% Prepare MUMC ME-GRE data for QSM and MIMM pipeline
%
% Converts per-echo NIfTI files (from dcm2niix) into 4D magnitude and
% phase arrays expected by run_QSM.m and run_MIMM_MUMC.m.
%
% Philips phase scaling (from JSON sidecar):
%   Stored as INT16. niftiread applies scl_slope/scl_inter automatically,
%   returning values in milli-radians. Divide by 1000 to get radians.
%
% Input file convention:
%   <prefix>-ME_GRE_e1.nii.gz       echo 1 magnitude
%   <prefix>-ME_GRE_e1_ph.nii.gz    echo 1 phase
%   ...
%   <prefix>-ME_GRE_e5.nii.gz       echo 5 magnitude
%   <prefix>-ME_GRE_e5_ph.nii.gz    echo 5 phase
%
% Outputs (written to subj_dir):
%   magnitude.nii.gz    [224 x 224 x 154 x 5]  double, RWV-scaled
%   phase.nii.gz        [224 x 224 x 154 x 5]  double, radians [-pi, pi]

%% --- Paths (adapt per subject) ---

data_dir  = '/home/tijn-saes/Documents/Internship/ME_GRE/';
subj_dir  = '/home/tijn-saes/Documents/Internship/ME_GRE/';
prefix    = '501';   % subject/series number prefix

n_echoes  = 5;
TE        = [0.006001, 0.012, 0.018, 0.024, 0.030];   % seconds, from JSON

%% --- Load and stack echoes ---

mag_4d = [];
pha_4d = [];

for e = 1:n_echoes
    mag_file      = fullfile(data_dir, sprintf('%s-ME_GRE_e%d.nii.gz',       prefix, e));
    pha_file      = fullfile(data_dir, sprintf('%s-ME_GRE_e%d_ph.nii.gz',    prefix, e));
    mag_json_file = fullfile(data_dir, sprintf('%s-ME_GRE_e%d.json',         prefix, e));
    pha_json_file = fullfile(data_dir, sprintf('%s-ME_GRE_e%d_ph.json',      prefix, e));

    % Read raw integer values (dcm2niix does not embed Philips scaling in NIfTI header)
    mag_vol = double(niftiread(mag_file));
    pha_vol = double(niftiread(pha_file));

    % Apply Philips RWV scaling from JSON sidecar
    mag_json = jsondecode(fileread(mag_json_file));
    pha_json = jsondecode(fileread(pha_json_file));

    mag_vol = mag_vol * mag_json.PhilipsRWVSlope + mag_json.PhilipsRWVIntercept;

    % Phase: RWV scaling gives milli-radians → divide by 1000 for radians
    pha_vol = (pha_vol * pha_json.PhilipsRWVSlope + pha_json.PhilipsRWVIntercept) / 1000;

    if isempty(mag_4d)
        sz     = size(mag_vol);
        mag_4d = zeros([sz, n_echoes]);
        pha_4d = zeros([sz, n_echoes]);
    end

    mag_4d(:,:,:,e) = mag_vol;
    pha_4d(:,:,:,e) = pha_vol;

    fprintf('Echo %d loaded: mag [%.1f %.1f]  phase [%.3f %.3f] rad\n', e, ...
        min(mag_vol(:)), max(mag_vol(:)), min(pha_vol(:)), max(pha_vol(:)));
end

%% --- Sanity check: phase should be in [-pi, pi] ---

if max(abs(pha_4d(:))) > pi * 1.01
    warning('Phase values exceed [-pi, pi] after scaling (max |val| = %.3f rad). Check Philips scaling.', ...
        max(abs(pha_4d(:))));
else
    disp('Phase range OK: within [-pi, pi] rad.');
end

%% --- Save 4D NIfTIs ---

ref_info = niftiinfo(fullfile(data_dir, sprintf('%s-ME_GRE_e1.nii.gz', prefix)));

% Magnitude
mag_info                  = ref_info;
mag_info.Datatype         = 'double';
mag_info.BitsPerPixel     = 64;
mag_info.ImageSize        = size(mag_4d);
mag_info.PixelDimensions  = [ref_info.PixelDimensions(1:3), TE(2) - TE(1)];
niftiwrite(mag_4d, fullfile(subj_dir, 'magnitude.nii.gz'), mag_info, 'Compressed', true);
disp('Written: magnitude.nii.gz');

% Phase
pha_info                  = mag_info;
niftiwrite(pha_4d, fullfile(subj_dir, 'phase.nii.gz'), pha_info, 'Compressed', true);
disp('Written: phase.nii.gz');

fprintf('Done. Data shape: %s\n', mat2str(size(mag_4d)));
fprintf('Echo times: %s ms\n', mat2str(TE * 1000));
