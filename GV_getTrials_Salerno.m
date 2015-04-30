function [Trials] = GV_getTrials_Salerno(Data,Labels,params)
% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl

Trials(size(Labels,1))      = struct('Data', [],'Condition',[],'Run',[],'HRF',[],...
    'Onset',[], 'TimeWindow',[],'TimeWindowBaseline',[]);
Intervals                   = params.data.Intervals;

for ind = 1:size(Labels,1)
    ThisTrial               = Intervals(ind,:);
    Trials(ind).Condition   = ThisTrial(3);
    Trials(ind).Run         = ThisTrial(4);
    Trials(ind).Onset       = ThisTrial(1);
    TimeWindow              = ThisTrial(1) + (-params.fe.preFit:params.fe.postFit);
    TimeWindowBaseline      = ThisTrial(1) + (-params.fe.preFit:0);
    TrialData               = Data(TimeWindow,:);
    BaselineData            = Data(TimeWindowBaseline,:);
    switch lower(params.trial.normalization)
        case 'percsig'
           TrialData        = (TrialData-repmat(mean(BaselineData),size(TrialData,1),1))./repmat(mean(BaselineData),size(TrialData,1),1);
        case 'diff'
           TrialData        = (TrialData-repmat(mean(BaselineData),size(TrialData,1),1)); 
    end
    Trials(ind).Data        = TrialData;
    Trials(ind).HRF         = params.fe.StdHrf;
    Trials(ind).TimeWindow  = TimeWindow;
    Trials(ind).TimeWindowBaseline = TimeWindowBaseline;
end