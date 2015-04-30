function [Ealpha, Ebeta, Omega, b, Evarphi, tau, nu] = GR_ProvaScriptGirolami(KernelBasis, KernelPoints,T)
%GR_PROVASCRIPTGIROLAMI Summary of this function goes here
%   Detailed explanation goes here
dim=size(KernelPoints(1).Kernel);
num=size(KernelPoints,1);
Ksupplied=zeros(dim(1),dim(2),num);
for i=1:8
    Ksupplied(:,:,i)=KernelPoints(i).Kernel;
end;

    %Ksupplied(:,:,1)=KernelBasis(1).stmp;
    %Ksupplied(:,:,2)=KernelBasis(1).dstmp;
    %Ksupplied(:,:,3)=KernelBasis(1).ddstmp;

%fprintf('\nSetting up preferences...');
%Set-up hyper-parameters
hyper.sigma = 3; %3     shape
hyper.varsigma = 50; %20 scale   molto piccata
hyper.tau = 3; %3       shape
hyper.nu = 50;   %20     scale
hyper.update = 1;
MAXIT = 100;
KWTol = 1e-4;
SampPars.NoS = 1000;
SampPars.MAXIT = 2;
SampPars.TOL = 1e-4;

Kdef2(1).ktype = 'supplied';%Because kernel is supplied
bias = 0;
display = 2;%Just iteration text

figure(2);
fprintf('\nRunning Algorithm...');
[Ealpha, Ebeta, Omega, b, Evarphi, tau, nu] = vb_classification_dirichlet(0,T,hyper,Kdef2,bias,display,MAXIT,KWTol,SampPars,Ksupplied);

end

