function [DataKernel] = GR_extractFeatures3KernelData_Salerno( Trials,HRFBasis )
%GV_EXTRACTFEATURESKERNELDATA_SALERNO Summary of this function goes here
%   Detailed explanation goes here

    
    dim=numel(Trials);
    numvoxel=size(Trials(1).Data,2);

    TempData1=zeros(dim,numvoxel);
    TempData2=zeros(dim,numvoxel);
    TempData3=zeros(dim,numvoxel);
    
    DataKernel = repmat(struct('stmp',zeros(dim,numvoxel),'dstmp',zeros(dim,numvoxel),'ddstmp',zeros(dim,numvoxel)), 1, 1);

    parfor i=1:dim

        TempData1(i,:)=GV_fitSingleTrial_Salerno(Trials(i).Data,(HRFBasis.stmp)');
        TempData2(i,:)=GV_fitSingleTrial_Salerno(Trials(i).Data,(HRFBasis.dstmp)');
        TempData3(i,:)=GV_fitSingleTrial_Salerno(Trials(i).Data,(HRFBasis.ddstmp)');
    end;
    
    DataKernel.stmp=TempData1;
    DataKernel.dstmp=TempData2;
    DataKernel.ddstmp=TempData3;
    
    
end

