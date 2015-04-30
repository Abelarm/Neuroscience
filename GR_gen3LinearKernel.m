function [ KernelBasis ] = GR_gen3LinearKernel( DataKernel )
%GR_GEN3LINEARKERNEL Summary of this function goes here
%   Detailed explanation goes here

    dim=size(DataKernel.stmp,1);
    
    KernelBasis = repmat(struct('stmp',zeros(dim,dim),'dstmp',zeros(dim,dim+1),'ddstmp',zeros(dim,dim+1)), 1, 1);

    Temp= (DataKernel.stmp)*(DataKernel.stmp)';
    q=1./(diag(Temp));
    Temp=[Temp .* repmat(q',size(DataKernel(1).stmp,1),1) repmat(1,size(DataKernel(1).stmp,1),1)];
    
    KernelBasis.stmp=Temp;
    
    Temp= (DataKernel.dstmp)*(DataKernel.dstmp)';
    q=1./(diag(Temp));
    Temp=[Temp .* repmat(q',size(DataKernel(1).dstmp,1),1) repmat(1,size(DataKernel(1).dstmp,1),1)];
    
    KernelBasis.dstmp=Temp;
    
    Temp= (DataKernel.ddstmp)*(DataKernel.ddstmp)';
    q=1./(diag(Temp));
    Temp=[Temp .* repmat(q',size(DataKernel(1).ddstmp,1),1) repmat(1,size(DataKernel(1).ddstmp,1),1)];
    
    KernelBasis.ddstmp=Temp;
    
    
end

