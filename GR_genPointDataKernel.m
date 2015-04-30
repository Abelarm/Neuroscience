function [ DataKernel ] = GR_genPointDataKernel( Trials )
%GENPOINTKERNEL Summary of this function goes here
%   Detailed explanation goes here



    sizesTrials=size(Trials,2);

    TempKernel= zeros(numel(Trials),size(Trials(1).Data,2));
    points=size(Trials(1).Data,1);
    

    DataKernel = repmat(struct('Data',zeros(numel(Trials),size(Trials(1).Data,2))), points, 1);

    for k=1:points
        DataKernel(k)=struct('Data',zeros(numel(Trials),size(Trials(1).Data,2)));
    end;

    for k=1:points
        for i=1:sizesTrials
            TempKernel(i,:)=Trials(i).Data(k,:);
        end;
        DataKernel(k).Data=TempKernel;
    end;


end

