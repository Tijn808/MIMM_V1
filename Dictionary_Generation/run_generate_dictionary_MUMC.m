% Generate MIMM dictionary for MUMC ME-GRE acquisition protocol.
% TEs: [6, 12, 18, 24, 30] ms — no interpolation needed in MIMM matching.
%
% Runtime: ~38 hours (same as original paper)
% Output:  Dictionary/MIMM_dictionary_stochastic_MUMC.mat
%
% Run with:
%   matlab -batch "run('/path/to/run_generate_dictionary_MUMC.m')"

mimm_root = fileparts(fileparts(mfilename('fullpath')));

addpath(fullfile(mimm_root, 'Dictionary_Generation'));

out_file = fullfile(mimm_root, 'Dictionary', 'MIMM_dictionary_stochastic_MUMC.mat');

fprintf('Starting MUMC dictionary generation...\n');
fprintf('TEs: [6 12 18 24 30] ms\n');
fprintf('Output: %s\n\n', out_file);

tic
dictionary = generate_dictionary_stochastic_MUMC();
elapsed = toc;

fprintf('\nDone. Elapsed time: %.1f hours\n', elapsed/3600);
fprintf('Dictionary size: %d entries\n', size(dictionary.signal, 1));

save(out_file, 'dictionary', '-v7.3');
fprintf('Saved to: %s\n', out_file);
