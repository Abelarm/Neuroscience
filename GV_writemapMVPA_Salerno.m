function [NameSave] = GV_writemapMVPA_Salerno(params,namefiles,results,MaskSize)
% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl

NameSave                        = [namefiles.NameSave '_' params.fe.FeatExt];


% NameMask                        = namefiles.NameMask;
% msk                             = xff(NameMask);
% Resolution                      = msk.Resolution;
% Mask                            = msk.Mask;
% MaskSize                        = size(msk.Mask);
% msk.ClearObject;



cmap                            = cell(1);
map                             = zeros(MaskSize);
map(params.data.vox)            = zscore(mean(results.temp_model_weights,2));
cmap{1}                         = struct('LowerThreshold',1.2,'UpperThreshold',10,'Name','Predictive map','CMPData',map,...
    'ParamValues', 1-results.mean_err);



ic                              = xff('ica');
ic.Resolution                   = 2;
ic.NrOfMapParameters            = 1;
ic.ShowParamsRangeFrom          = 0;
ic.ShowParamsRangeTo            = 1;
ic.Map(1).CMPData               = cmap{1}.CMPData;
ic.Map(1).LowerThreshold        = cmap{1}.LowerThreshold;
ic.Map(1).UpperThreshold        = cmap{1}.UpperThreshold;
ic.Map(1).Name                  = cmap{1}.Name;
ic.Map(1).EnableClusterCheck    = 1;
ic.Map(1).UseRGBColor           = 0;
ic.MapParameter(1).Name         = 'Accuracy';
ic.MapParameter(1).Values       = 1-results.mean_err;

ic.SaveAs([NameSave '.ica']);
ic.ClearObject;




