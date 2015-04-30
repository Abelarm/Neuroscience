function [X,params] = importfMRIdata(NameFiles,params)
% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl


NameVtc             = NameFiles.NameVtc;
NameMask            = NameFiles.NameMask;



% reading the mask and selecting the voxels within the mask
msk                 = xff(NameMask);
VoxelsUsed          = find(msk.Mask(:)>0);
msk.ClearObject;



% reading data from vtc and protocol
X                   = [];   % data matrix
MeanRuns            = zeros(numel(NameVtc),numel(VoxelsUsed));
StdRuns             = zeros(numel(NameVtc),numel(VoxelsUsed));
NumVol       = zeros(numel(NameVtc),1);

for idx             = 1:numel(NameVtc)
   
    v               = xff(NameVtc{idx});
    if idx == 1
        params.data.TR          = v.TR;
    end
    % reading the data
    ThisData        = double(v.VTCData(:,VoxelsUsed));
    NumVol(idx)     = v.NrOfVolumes;
    v.ClearObject;
    % storing mean and standard deviation of each run
    MeanRuns(idx,:) = mean(ThisData,1);
    StdRuns(idx,:)  = std(ThisData,0,1);
    
    % remove mean of each voxel if asked
    if params.pre_process.rem_meanvtc
        ThisData    = bsxfun(@minus,ThisData,MeanRuns(idx,:));
    end
    % remove standard deviation of each voxel if asked
    if params.pre_process.rem_varvtc
        ThisData    = bsxfun(@rdivide,ThisData,StdRuns(idx,:));
    end
    
    X = [X; ThisData];
end


IndexSelMean        = find(sum(MeanRuns > params.pre_process.mean_threshold,1) == numel(NameVtc));
IndexSelStd         = find(sum(StdRuns  < params.pre_process.mean_threshold,1) == numel(NameVtc));

IndexSel            = intersect(IndexSelMean,IndexSelStd);

X                   = X(:,IndexSel);
VoxelsUsed          = VoxelsUsed(IndexSel);
params.data.vox     = VoxelsUsed;
params.data.NumVol  = NumVol;





