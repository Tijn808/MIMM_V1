mimm_root = fileparts(fileparts(mfilename('fullpath')));

orig = load(fullfile(mimm_root, 'Dictionary', 'MIMM_dictionary_stochastic.mat'), 'dictionary');
mumc = load(fullfile(mimm_root, 'Dictionary', 'MIMM_dictionary_stochastic_MUMC.mat'), 'dictionary');

fprintf('Original TEs (ms): %s\n', mat2str(orig.dictionary.TE * 1e3));
fprintf('MUMC TEs (ms):     %s\n', mat2str(mumc.dictionary.TE * 1e3));
fprintf('Original signal size: %s\n', mat2str(size(orig.dictionary.signal)));
fprintf('MUMC signal size:     %s\n', mat2str(size(mumc.dictionary.signal)));

% Check values match exactly at the subsampled indices
te_idx = [3 5 7 9 11];
orig_sub = orig.dictionary.signal(:, te_idx);
mumc_sig = mumc.dictionary.signal;
max_diff = max(abs(orig_sub(:) - mumc_sig(:)));
fprintf('Max difference between subsampled original and MUMC dict: %.2e\n', max_diff);

if max_diff < 1e-10
    disp('PASS: MUMC dictionary is an exact subsample of the original.');
else
    disp('FAIL: Values do not match exactly — something went wrong.');
end
