function [ TrialsMed ] = GR_medTrials( Trials )
%GR_MEDTRIALS Summary of this function goes here
%   Detailed explanation goes here

  TrialsMed=Trials;
  num=numel(TrialsMed);
  dim=size(TrialsMed(1).Data);
  
    parfor i=1:num
        for j=1:dim(2)
            media=mean(Trials(i).Data(:,j))
            TrialsMed(i).Data(:,j)=Trials(i).Data(:,j)-media;
        end;
    end;

end

