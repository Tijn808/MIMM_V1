% Prepare MUMC ME-GRE data for QSM and MIMM pipeline
%
% Converts per-echo NIfTI files (from dcm2niix) into 4D magnitude and
% phase arrays expected by run_QSM_chisep.m and run_MIMM_MUMC.m.
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

%% --- Paths (loaded from MUMC_pipeline/paths.m) ---

if ~exist('mimm_root', 'var')
    paths_file = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'paths.m');
    if ~exist(paths_file, 'file')
        error('paths.m not found. Copy MUMC_pipeline/paths_template.m to paths.m and fill in your paths.');
    end
    run(paths_file);
end
data_dir  = input_dir;
subj_dir  = input_dir;
if ~exist('prefix', 'var')
    prefix = '501';   % adapt per subject, or supply via run_subject.m / run_cohort.sh
end

%% --- MUMC format (gremag.nii / grepha.nii: pre-stacked 4D, from MUMC 010) ---
% MUMC's 010_DicomToNifti writes the ME-GRE as two 4D files (all echoes in one),
% uncompressed, instead of five separate <prefix>-ME_GRE_e<n>.nii.gz echoes.
% If we see that layout, read it directly and skip the per-echo stacking below.
% Look both in data_dir and in data_dir/nifti (010 writes into a nifti/ subdir).
mumc_mag = '';
for dd = {data_dir, fullfile(data_dir, 'nifti')}
    for c = {'gremag.nii', 'gremag.nii.gz'}
        f = fullfile(dd{1}, c{1});
        if exist(f, 'file'); mumc_mag = f; break; end
    end
    if ~isempty(mumc_mag); break; end
end
if ~isempty(mumc_mag)
    mumc_pha = '';
    for dd = {data_dir, fullfile(data_dir, 'nifti')}
        for c = {'grepha.nii', 'grepha.nii.gz'}
            f = fullfile(dd{1}, c{1});
            if exist(f, 'file'); mumc_pha = f; break; end
        end
        if ~isempty(mumc_pha); break; end
    end
    if isempty(mumc_pha)
        error('Found %s but no grepha (phase) file in %s.', mumc_mag, data_dir);
    end
    fprintf('MUMC format detected: %s + %s\n', mumc_mag, mumc_pha);

    mag_4d = double(niftiread(mumc_mag));
    pha_4d = double(niftiread(mumc_pha));

    % Apply Philips RWV scaling from the JSON sidecars (same convention as the
    % per-echo path). One slope per 4D file; for Philips it is constant across
    % echoes. (MIMM matches the decay *shape*, so a global magnitude scale is
    % harmless; the phase must end up in radians.)
    for pair = { {mumc_mag, 'mag'}, {mumc_pha, 'pha'} }
        f    = pair{1}{1};
        jf   = regexprep(f, '\.nii(\.gz)?$', '.json');
        if exist(jf, 'file')
            j = jsondecode(fileread(jf));
            if isfield(j, 'PhilipsRWVSlope')
                s = j.PhilipsRWVSlope; b = 0;
                if isfield(j, 'PhilipsRWVIntercept'); b = j.PhilipsRWVIntercept; end
                if strcmp(pair{1}{2}, 'mag'); mag_4d = mag_4d*s + b; else; pha_4d = pha_4d*s + b; end
            end
        end
    end
    pha_4d = pha_4d / 1000;   % milli-radians -> radians (Philips convention)

    if max(abs(pha_4d(:))) > pi*1.01
        warning(['Phase outside [-pi,pi] (max |val|=%.3f). Rescaling to [-pi,pi]. ' ...
                 'Check the Philips phase scaling for this subject.'], max(abs(pha_4d(:))));
        pha_4d = pha_4d / max(abs(pha_4d(:))) * pi;
    else
        disp('Phase range OK: within [-pi, pi] rad.');
    end

    TE = [6.001 12 18 24 30];   % protocol TEs (ms); QSM/MIMM use their own TE too

    ref = niftiinfo(mumc_mag);
    ref.Datatype = 'double'; ref.BitsPerPixel = 64;
    ref.ImageSize = size(mag_4d);
    pd = ref.PixelDimensions;
    if numel(pd) < 4; pd = [pd(1:3), 6]; end
    ref.PixelDimensions = [pd(1:3), TE(2)-TE(1)];
    niftiwrite(mag_4d, fullfile(subj_dir, 'magnitude.nii.gz'), ref, 'Compressed', true);
    niftiwrite(pha_4d, fullfile(subj_dir, 'phase.nii.gz'),     ref, 'Compressed', true);
    fprintf('MUMC prepare done. Shape: %s, phase [%.3f, %.3f] rad\n', ...
        mat2str(size(mag_4d)), min(pha_4d(:)), max(pha_4d(:)));
    return
end

n_echoes  = 5;
TE        = zeros(1, n_echoes);   % filled from JSON below

%% --- Load and stack echoes ---

% Preallocate using echo-1 header so the analyser does not warn about growing arrays.
info_e1 = niftiinfo(fullfile(data_dir, sprintf('%s-ME_GRE_e1.nii.gz', prefix)));
mag_4d  = zeros([info_e1.ImageSize, n_echoes]);
pha_4d  = zeros([info_e1.ImageSize, n_echoes]);
clear info_e1

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

    TE(e) = mag_json.EchoTime;   % seconds, as reported by scanner

    mag_vol = mag_vol * mag_json.PhilipsRWVSlope + mag_json.PhilipsRWVIntercept;

    % Phase: RWV scaling gives milli-radians → divide by 1000 for radians
    pha_vol = (pha_vol * pha_json.PhilipsRWVSlope + pha_json.PhilipsRWVIntercept) / 1000;

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
