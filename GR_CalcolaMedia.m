function [ medErr ] = GR_CalcolaMedia( ERR,dim2 )
%GR_CALCOLAMEDIA Summary of this function goes here
%   Detailed explanation goes here

    dim=size(ERR);
    medErr=0;
    
    for i=1:dim(2)
        if(mod(i,10)==10)
            medErr = medErr + ERR(i)*24;
        else
            medErr = medErr + ERR(i)*8;
        end;
    end;
    
    medErr=medErr/(dim2*5);


end

