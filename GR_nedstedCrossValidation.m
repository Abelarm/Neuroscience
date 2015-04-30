function [ ERRORITOTALI ] = GR_nedstedCrossValidation( Kernel,T,Splits )
%GR_NEDSTEDCROSSVALIDATION Summary of this function goes here
%   Detailed explanation goes here

MAT=repmat(struct('val',zeros(2)),4,4);
MAT(1,1).val=[1,0.5];
MAT(1,2).val=[1,1];
MAT(1,3).val=[1,3];
MAT(1,4).val=[1,5];
MAT(2,1).val=[10,0.5];
MAT(2,2).val=[10,1];
MAT(2,3).val=[10,3];
MAT(2,4).val=[10,5];
MAT(3,1).val=[30,0.5];
MAT(3,2).val=[30,1];
MAT(3,3).val=[30,3];
MAT(3,4).val=[30,5];
MAT(4,1).val=[50,0.5];
MAT(4,2).val=[50,1];
MAT(4,3).val=[50,3];
MAT(4,4).val=[50,5];
ERRORITOTALI=zeros(1,numel(Splits));


NITER=1;

parfor k=10:numel(Splits)
    SubSplit=GR_split85_73(Splits{k}.indTrain);
    ERR=zeros(16,size(SubSplit,2));
    t=1;
    for i=1:size(MAT,1)
        for j=1:size(MAT,2)
            %Param(1)=MAT(i,j).val(1);
            %Param(2)=MAT(i,j).val(2);
            fprintf('Valutando i parametri %f %f\n',MAT(i,j).val(1), MAT(i,j).val(2));
            for q=1:size(SubSplit,2)
                ERR(t,q)  = GR_classifyVal3Kernel(Kernel,SubSplit{q},T,[MAT(i,j).val(1) MAT(i,j).val(2)],NITER);
                
            end;
            t=t+1;
            
        end;
    end;
    TotErr=sum(ERR);
    [minimo,ind]=min(TotErr);
    [x,y]=ind2sub([4 4],ind);
    %Param(1)=MAT(x,y).val(1);
    %Param(2)=MAT(x,y).val(2);
    ERRORITOTALI(k)=GR_classifyVal3Kernel(Kernel,Splits{k},T,[MAT(i,j).val(1) MAT(i,j).val(2)],NITER);
    fprintf('Valutato %d° split\n',k);
end;

end

