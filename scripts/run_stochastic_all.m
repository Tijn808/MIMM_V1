cd(fileparts(fileparts(mfilename('fullpath'))))
MIMM_set_path
stoch = load('MIMM_dictionary_stochastic.mat');
dict = stoch.dictionary;
lambda_chi = 0.015;
load('FA.mat')
load('QSM.mat')
load('theta.mat')
load('iField.mat','iField','TE')

% Basic MIMM
orientation_strategy = 'basic';
MIMM_basic = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, orientation_strategy);
disp('Basic done!')

% DTI Orientation Informed MIMM
orientation_strategy = 'orientation_informed';
MIMM_DTI = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, orientation_strategy, FA_DTI, theta_DTI);
disp('DTI done!')

% Atlas Orientation Informed MIMM
orientation_strategy = 'orientation_informed';
MIMM_Atlas = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, orientation_strategy, FA_atlas, theta_atlas);
disp('Atlas done!')

save('Example_Results/stochastic_MIMM_results_new.mat', 'MIMM_basic', 'MIMM_DTI', 'MIMM_Atlas')
disp('All done!')
