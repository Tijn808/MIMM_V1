% Test pipeline on phantom data (no SEPIA required)
% Creates a threshold-based mask and dummy QSM, then runs MIMM Basic.
% Purpose: verify the full preprocessing → MIMM chain works on MUMC data format.

cd('/home/tijn-saes/Documents/Internship/MIMM')
MIMM_set_path

data_dir   = '/home/tijn-saes/Documents/Internship/ME_GRE/';
output_dir = '/home/tijn-saes/Documents/Internship/ME_GRE/test_output/';
if ~exist(output_dir, 'dir'); mkdir(output_dir); end

TE         = [0.006001, 0.012, 0.018, 0.024, 0.030];
lambda_chi = 0.015;

%% --- Load 4D magnitude ---
disp('Loading magnitude...')
mag_info = niftiinfo(fullfile(data_dir, 'magnitude.nii.gz'));
iField   = double(niftiread(mag_info));   % [224 x 224 x 154 x 5]

%% --- Create threshold mask (replaces BET for phantom) ---
disp('Creating mask...')
mag_e1     = iField(:,:,:,1);
threshold  = 0.1 * max(mag_e1(:));
Brain_Mask = mag_e1 > threshold;

% Save mask
ref        = mag_info;
ref.ImageSize      = size(Brain_Mask);
ref.Datatype       = 'uint8';
ref.BitsPerPixel   = 8;
ref.PixelDimensions = mag_info.PixelDimensions(1:3);
niftiwrite(uint8(Brain_Mask), fullfile(output_dir, 'mask.nii.gz'), ref, 'Compressed', true);
fprintf('Mask: %d voxels active (%.1f%%)\n', sum(Brain_Mask(:)), 100*mean(Brain_Mask(:)));

%% --- Dummy QSM (zeros — no SEPIA available yet) ---
QSM = zeros(size(Brain_Mask));

%% --- Load dictionary ---
disp('Loading dictionary...')
stoch = load('Dictionary/MIMM_dictionary_stochastic.mat');
dict  = stoch.dictionary;

%% --- Run MIMM Basic ---
disp('Running MIMM Basic...')
tic
MIMM_basic = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, 'basic');
t = toc;
fprintf('MIMM done in %.1f seconds.\n', t)

%% --- Save MVF as NIfTI ---
mvf_info              = ref;
mvf_info.Datatype     = 'double';
mvf_info.BitsPerPixel = 64;
niftiwrite(MIMM_basic.MVF, fullfile(output_dir, 'MVF_basic.nii.gz'), mvf_info, 'Compressed', true);

fprintf('MVF range: [%.4f, %.4f]\n', min(MIMM_basic.MVF(:)), max(MIMM_basic.MVF(:)))
disp('Done. Output written to test_output/')
