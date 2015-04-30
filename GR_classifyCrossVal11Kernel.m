function [ err ] = GR_classifyCrossVal11Kernel( Kernel,KernelBasis,Splits, T,Param,NITER )
%GR_CLASSIFYCROSSVAL_SALERNO Summary of this function goes here
%   Detailed explanation goes here
    
    err=0;
    nsplit= numel(Splits);
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
    
    for i=1:nsplit
        TrainKernel=zeros(size(Splits{i}.indTrain,1),size(Splits{i}.indTrain,1),dim(1)+3);
        TestKernel=zeros(size(Splits{i}.indTest,1),size(Splits{i}.indTrain,1),dim(1)+3);
            parfor k=1:dim(1)
                Kernel(k).Kernel(Splits{i}.indTrain,Splits{i}.indTrain);
                TrainKernel(:,:,k)=Kernel(k).Kernel(Splits{i}.indTrain,Splits{i}.indTrain);
                TestKernel(:,:,k)=Kernel(k).Kernel(Splits{i}.indTest,Splits{i}.indTrain);
            end;  
            
                KernelBasis(1).stmp(Splits{i}.indTrain,Splits{i}.indTrain);
                TrainKernel(:,:,9)=KernelBasis(1).stmp(Splits{i}.indTrain,Splits{i}.indTrain);
                TestKernel(:,:,9)=KernelBasis(1).stmp(Splits{i}.indTest,Splits{i}.indTrain);
                
                KernelBasis(1).dstmp(Splits{i}.indTrain,Splits{i}.indTrain);
                TrainKernel(:,:,10)=KernelBasis(1).dstmp(Splits{i}.indTrain,Splits{i}.indTrain);
                TestKernel(:,:,10)=KernelBasis(1).dstmp(Splits{i}.indTest,Splits{i}.indTrain);
                
                KernelBasis(1).dstmp(Splits{i}.indTrain,Splits{i}.indTrain);
                TrainKernel(:,:,11)=KernelBasis(1).ddstmp(Splits{i}.indTrain,Splits{i}.indTrain);
                TestKernel(:,:,11)=KernelBasis(1).ddstmp(Splits{i}.indTest,Splits{i}.indTrain);
            
            [Ealpha, Ebeta, Omega, b, Evarphi, tau, nu] = vb_classification_dirichlet(0,T(Splits{i}.indTrain),hyper,Kdef2,bias,display,MAXIT,KWTol,SampPars,TrainKernel);
            
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
                    
                valReal=T(Splits{i}.indTest(q));
                
                if(valPred~=valReal)
                    errPar=errPar+1;
                end;
            end;
            errPar=errPar/size(Splits{i}.indTest,1);
            err(i)=errPar;
        end;

end

