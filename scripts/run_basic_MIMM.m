cd('/home/tijn-saes/Documents/Internship/MIMM')
MIMM_set_path
stoch = load('MIMM_dictionary_stochastic.mat');
dict = stoch.dictionary;
lambda_chi = 0.015;
load('FA.mat')
load('QSM.mat')
load('theta.mat')
load('iField.mat','iField','TE')
orientation_strategy = 'basic';
MIMM_basic = MIMM(dict, lambda_chi, QSM, Brain_Mask, iField, TE, orientation_strategy);
save('Example_Results/stochastic_MIMM_basic.mat','MIMM_basic')
disp('Done!')
