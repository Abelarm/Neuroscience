function [HRFBasis] = GR_extractHRFBasis_Salerno( )
%GV_EXTRACTHRFBASIS_SALERNO Summary of this function goes here
%   Detailed explanation goes here

    %system('python BasisEstimations.py');
    load('basis.mat')
    
    HRFBasis=struct('stmp',zeros(8,1));
    HRFBasis=struct('dstmp',zeros(8,1));
    HRFBasis=struct('ddstmp',zeros(8,1));
    
    bas1med=mean(basis(1,:));
    bas2med=mean(basis(2,:));
    bas3med=mean(basis(3,:));
    HRFBasis.stmp=(basis(1,:)-bas1med);
    HRFBasis.dstmp=(basis(2,:)-bas2med);
    HRFBasis.ddstmp=(basis(3,:)-bas3med);
    
    %delete('basis.mat')

end

