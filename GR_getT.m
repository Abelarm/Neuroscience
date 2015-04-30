function [ T ] = GR_getT( Trials )
%GR_GETT Summary of this function goes here
%   Detailed explanation goes here

    dim=size(Trials,2);
    T=zeros(dim,1);
    
    for i=1:dim
        if Trials(i).Condition==2
            T(i)=1;
        else
            T(i)=-1;
        end;
    end;
end

