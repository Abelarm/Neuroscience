function [DM,Labels,params]      = GV_importProtocol_Salerno(NameFiles,params)
% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl

NamePrt             = NameFiles.NamePrt;
hrfval              = params.hrf;


DM                  = [];     % Design Matrix
Intervals           = [];


for idx             = 1:numel(NamePrt)
    
    p = xff(NamePrt{idx});

    
    ThisDesign      = p.CreateSDM(struct('hshape','twogamma','hpttp',hrfval.hpttp,'hnttp',hrfval.hnttp,'hpnr',hrfval.hpnr,...
        'hons',hrfval.hons,'hpdsp',hrfval.hpdsp,'hndsp',hrfval.hndsp,'prtr',...
        params.data.TR , 'nvol', params.data.NumVol(idx), 'rcond', 1));
    DM              = [DM;[ThisDesign.RTCMatrix  ones(params.data.NumVol(idx,1),1)*idx ]];
    
    
    % creating the "intervals" array
    p.ConvertToVol(params.data.TR);
    switch lower(params.Design)
        case {'block','slow'}
            
            OnOffsets           = [];
            for idy             = 1:p.NrOfConditions
                z               = p.Cond(idy).OnOffsets;
                z(z>params.data.NumVol(idx)) = params.data.NumVol(idx);
                z               = z + sum(params.data.NumVol(1:idx-1));
                if ~isempty(z)
                    OnOffsets   = [OnOffsets ; [z  idy*ones(size(z,1),1) idx*ones(size(z,1),1)]];
                end
            end
            Intervals           = [Intervals; OnOffsets];
    end

    
end
% removing the first condition, it is rest! (only for this experiment...)
Intervals(Intervals(:,3)==1,:)  = [];
params.data.Intervals           = Intervals;
Labels                = Intervals(:,[3 4]);

% getting the standard hrf for the given window
[~,FirstCond]                   =  min(Intervals(:,1));
HrfTrial                        = DM(Intervals(FirstCond,1) + (-params.fe.preFit:params.fe.postFit),Intervals(FirstCond,3)-1);
params.fe.StdHrf                = HrfTrial;
