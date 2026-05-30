% Subsample original MIMM dictionary to MUMC echo times.
% MUMC TEs [6,12,18,24,30] ms are exact grid points in the original
% dictionary (0:3:60 ms), so no interpolation is needed — just extract
% the matching columns directly.

mimm_root = fileparts(fileparts(mfilename('fullpath')));

orig_file = fullfile(mimm_root, 'Dictionary', 'MIMM_dictionary_stochastic.mat');
out_file  = fullfile(mimm_root, 'Dictionary', 'MIMM_dictionary_stochastic_MUMC.mat');

stoch = load(orig_file, 'dictionary');
dictionary = stoch.dictionary;

% Original TEs: 0:3:60 ms (21 points)
TE_orig = dictionary.TE * 1e3;   % convert to ms
TE_mumc = [6 12 18 24 30];       % MUMC echo times in ms

% Find indices
idx = arrayfun(@(t) find(abs(TE_orig - t) < 1e-6), TE_mumc);
fprintf('Subsampling at indices: %s\n', mat2str(idx));
fprintf('Corresponding TEs: %s ms\n', mat2str(TE_orig(idx)));

% Subsample signal matrix
dictionary.signal = dictionary.signal(:, idx);
dictionary.TE     = TE_mumc * 1e-3;   % back to seconds

fprintf('Original signal size: 20000 x 21\n');
fprintf('Subsampled signal size: %s\n', mat2str(size(dictionary.signal)));

save(out_file, 'dictionary', '-v7.3');
fprintf('Saved to: %s\n', out_file);
