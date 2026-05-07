% Export MIMM example results to NIfTI
%
% The example data (Zenodo) is stored as .mat files, so there is no
% existing NIfTI header to borrow. A minimal header is constructed from
% the known acquisition parameters (1 mm isotropic, 224x224x154).
%
% Outputs written to: Example_Results/nifti/

cd('/home/tijn-saes/Documents/Internship/MIMM')
MIMM_set_path

%% --- Load results ---

disp('Loading stochastic results...')
load('Example_Results/stochastic_MIMM_results.mat', 'MIMM_basic', 'MIMM_DTI', 'MIMM_Atlas')

%% --- Build reference NIfTI header ---
% niftiinfo() requires a file, so bootstrap by writing a dummy volume first.

voxel_mm = [1, 1, 1];

dummy_file = fullfile(tempdir, 'mimm_ref.nii');
niftiwrite(single(MIMM_basic.MVF), dummy_file);
ref = niftiinfo(dummy_file);
delete(dummy_file);

ref.PixelDimensions  = voxel_mm;
ref.Datatype         = 'double';
ref.BitsPerPixel     = 64;
ref.SpaceUnits       = 'Millimeter';

%% --- Output directory ---

out_dir = 'Example_Results/nifti';
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

%% --- Export all fields for each variant ---

variants      = {MIMM_basic,  MIMM_DTI,  MIMM_Atlas};
variant_names = {'basic',     'DTI',     'Atlas'};

fields = {'MVF', 'FVF', 'g_ratio', 'R2s', 'chi_iron_est', 'chi_myelin', ...
          'theta_est', 'error'};

for vi = 1:numel(variants)
    mv    = variants{vi};
    vname = variant_names{vi};

    sub_dir = fullfile(out_dir, vname);
    if ~exist(sub_dir, 'dir'); mkdir(sub_dir); end

    for fi = 1:numel(fields)
        fn   = fields{fi};
        data = double(mv.(fn));

        info             = ref;
        info.ImageSize   = size(data);

        fname = fullfile(sub_dir, [fn '.nii.gz']);
        niftiwrite(data, fname, info, 'Compressed', true);
        fprintf('  Written: %s/%s/%s.nii.gz\n', out_dir, vname, fn);
    end
end

disp('Export complete.')
fprintf('Files written to: %s\n', fullfile(pwd, out_dir))
