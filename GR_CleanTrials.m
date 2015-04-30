function [ Trial2,Label2 ] = GR_CleanTrials( TrialOld, LabelOld )
%NEW TRIAL Summary of this function goes here
%   Detailed explanation goes here
    sizeT=size(TrialOld,2);
    Trial2=struct('Data', [],'Condition',[],'Run',[],'HRF',[],...
    'Onset',[], 'TimeWindow',[],'TimeWindowBaseline',[]);
    cnt=1;
    for i=1:sizeT
        if TrialOld(i).Condition==2 || TrialOld(i).Condition==3
            Trial2(cnt)=TrialOld(i);
            Label2(cnt,:)=LabelOld(i,:);
            cnt=cnt+1;
        end
    end
    %cnt=1;
    %for i=1:sizeT
    %    if LabelOld(i,1)==2 || LabelOld(i,1)==3
           % Label2(cnt,:)=LabelOld(i,:);
           % cnt=cnt+1;
   % end
    
end

