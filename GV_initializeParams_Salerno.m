function params = GV_initializeParams_Salerno
% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl

params.pre_process.rem_meanvtc          = 0;  % removing the mean of each voxel per run: may be useful if no trial normalization is performed
params.pre_process.rem_varvtc           = 0;    % normalizing the variance of each voxel per run: may be useful if no trial normalization is performed  
params.pre_process.mean_threshold       = 200; % thesholds for mean of  voxels marked for removal
params.pre_process.std_threshold        = 100; % thesholds for std of voxels marked for removal


% hemodynamic paramters (these are the standard ones)
params.hrf.hpttp                        = 5;
params.hrf.hnttp                        = 15;
params.hrf.hpnr                         = 6;
params.hrf.hons                         = 0;
params.hrf.hpdsp                        = 1;
params.hrf.hndsp                        = 1;

% design type: it can be either block, slow or fast
params.Design                           = 'slow';

% choosing which conditions should be compared, in this case it is
% straightforward, 2 and 3
params.classes.C1                       = 2;
params.classes.C2                       = 3;

% cross-validation scheme,
% - loo: leave one out
% - lro: leave run out
% - kfold: k-fold cross validation, with additional parameters 
%       numfolds and repetitions
params.cv.scheme                        = 'kfold';'loo';'lro';
params.cv.numfolds                      = 10; % used only if cv.scheme == kfold
params.cv.niter                         = 5; % used only if cv.scheme == kfold


% trial normalization
params.trial.normalization              = 'percsig';


% feature estimation and extraction parameters
params.fe.FeatExt                       = 't'; %can be avg, beta or t
params.fe.UseEmpHrf                     = 0;    % if 1 derive hrf from train data, otherwise use standard
params.fe.preFit                        = 1;
params.fe.postFit                       = 6;
params.fe.AvgWindow                     = 2:5;

% params.fe.oneHRF                        = 1; % if ==1, one hrf for all the voxels
% params.fe.voxFixHRF                     = 100;
% params.fe.shiftHRF                      = 0;
% params.fe.shiftHRFvals                  = [-0.25 0.25 0.125];
% params.fe.expandHRF                     = 0;
% params.fe.expandHRFvals                 = [0.8 1.2 0.1];
% params.fe.useAllTrainingTrials          = 1;


% % permutation test
% params.perm.doperm                      = 1;
% params.perm.nperm                       = 2e3;


