function  [X]   = GV_extractFeatures_Salerno(Trials,Splits,params)
% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl

% here is the function to modify! I included Splits as input parameter,
% since the fitting will be (slightly) different for each partition of the
% data in training and testing dataset.


% This is the default option: 
% fitting each trial with the given hrf, no need to repeat the fitting
% for different splits
X                   = zeros(numel(Trials),size(Trials(1).Data,2));

% if you create a different model per training dataset, then create a 3-D
% array where each page represents a split
% X                   = zeros(numel(Trials),size(Trials(1).Data,2),numel(Splits));


parfor indTrial        = 1:numel(Trials)
    ThisTrial       = Trials(indTrial);
    switch lower(params.fe.FeatExt)
        case 'avg'
            X(indTrial,:)       = mean(ThisTrial.Data(params.fe.AvgWindow,:),1);
        case 'beta'
            [X(indTrial,:)]     = GV_fitSingleTrial_Salerno(ThisTrial.Data,ThisTrial.HRF);
        case 't'
            [~,X(indTrial,:)]   = GV_fitSingleTrial_Salerno(ThisTrial.Data,ThisTrial.HRF);
    end
%     X(indTrial,:) = z;
end
