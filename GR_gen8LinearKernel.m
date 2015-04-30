function [ Kernels ] = GR_gen8LinearKernel( DataKernel )
%GV_GENLINEARKERNEL Summary of this function goes here
%   Detailed explanation goes here

    Kernels = repmat(struct('Kernel',zeros(size(DataKernel(1).Data,1),size(DataKernel(1).Data,1))), size(DataKernel,1), 1);

    for k=1:size(DataKernel,1)
        Kernels(k)=struct('Kernel',zeros(size(DataKernel(1).Data,1),size(DataKernel(1).Data,1)+1));
    end;

    dim=numel(DataKernel);

    for i=1:dim
        Kernel=(DataKernel(i).Data)*(DataKernel(i).Data)';
        
        q=1./(diag(Kernel));
        Kernel=[Kernel .* repmat(q',size(DataKernel(1).Data,1),1) repmat(1,size(DataKernel(1).Data,1),1)];
        Kernels(i).Kernel=Kernel;
    end;
end

