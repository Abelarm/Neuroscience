function [ err ] = GR_classifyVal3Kernel( Kernel,Split,T,Param,NITER )
%GR_CLASSIFYCROSSVAL_SALERNO Summary of this function goes here
%   Detailed explanation goes here

err=0;
dim=size(Kernel);


%Set-up hyper-parameters
hyper.sigma = Param(2); %3     shape
hyper.varsigma = Param(1); %20 scale   molto piccata
hyper.tau = Param(2); %3       shape
hyper.nu = Param(1);   %20     scale
hyper.update = 1;
MAXIT = NITER;
KWTol = 1e-4;
SampPars.NoS = 1000;
SampPars.MAXIT = 2;
SampPars.TOL = 1e-4;

Kdef2(1).ktype = 'supplied';%Because kernel is supplied
bias = 0;
display = 0;

TrainKernel(:,:,1)=Kernel(1).stmp(Split.indTrain,Split.indTrain);
TestKernel(:,:,1)=Kernel(1).stmp(Split.indTest,Split.indTrain);


TrainKernel(:,:,2)=Kernel(1).dstmp(Split.indTrain,Split.indTrain);
TestKernel(:,:,2)=Kernel(1).dstmp(Split.indTest,Split.indTrain);


TrainKernel(:,:,3)=Kernel(1).ddstmp(Split.indTrain,Split.indTrain);
TestKernel(:,:,3)=Kernel(1).ddstmp(Split.indTest,Split.indTrain);


[Ealpha, Ebeta, Omega, b, Evarphi, tau, nu] = vb_classification_dirichlet(0,T(Split.indTrain),hyper,Kdef2,bias,display,MAXIT,KWTol,SampPars,TrainKernel);

KtestVB = Ebeta(1).*TestKernel(:,:,1);
for k = 2:size(TestKernel,3)
    KtestVB = KtestVB + Ebeta(k).*TestKernel(:,:,k);
end;

%Pred=tanh(KtestVB*Ealpha);
Pred=(1./(1+exp(-KtestVB*Ealpha)));
errPar=0;

for q=1:size(Pred,1)
    if(Pred(q)>0.5)
        valPred=1;
    else
        valPred=-1;
    end;
    
    valReal=T(Split.indTest(q));
    
    if(valPred~=valReal)
        errPar=errPar+1;
    end;
end;
errPar=errPar/size(Split.indTest,1);
err=errPar;
%err=err+errPar;

%err=err/45;

end

