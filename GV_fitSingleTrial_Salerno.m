function [B,varargout] = GV_fitSingleTrial_Salerno(X,HRF)
% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl



% if additional predictors are given as covariates, put them from the
% second column of HRF on 


DM              = [HRF(:) ones(size(HRF,1),1)];
Betas           = DM\X;
B               = Betas(1,:);
if nargout >= 2
     CC             = inv(DM'*DM);
     VarErr         = var(X- DM*Betas);
     C              = zeros(1,size(HRF,2)+1);
     C(1)           = 1;
     T              = C*Betas./sqrt(VarErr*(C*CC*C'));
     varargout{1}   = T;
end

