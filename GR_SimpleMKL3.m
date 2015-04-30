function [ err ] = GR_SimpleMKL3( Kernel,T,Splits)
%GR_SIMPLEMKL Summary of this function goes here
%   Detailed explanation goes here




C = [10];
verbose=2;

options.algo='svmclass'; % Choice of algorithm in mklsvm can be either
% 'svmclass' or 'svmreg'
%------------------------------------------------------
% choosing the stopping criterion
%------------------------------------------------------
options.stopvariation=1; % use variation of weights for stopping criterion
options.stopKKT=0;       % set to 1 if you use KKTcondition for stopping criterion
options.stopdualitygap=0; % set to 1 for using duality gap for stopping criterion

%------------------------------------------------------
% choosing the stopping criterion value
%------------------------------------------------------
options.seuildiffsigma=1e-2;        % stopping criterion for weight variation
options.seuildiffconstraint=0.1;    % stopping criterion for KKT
options.seuildualitygap=0.01;       % stopping criterion for duality gap

%------------------------------------------------------
% Setting some numerical parameters
%------------------------------------------------------
options.goldensearch_deltmax=1e-1; % initial precision of golden section search
options.numericalprecision=1e-8;   % numerical precision weights below this value
% are set to zero
options.lambdareg = 1e-8;          % ridge added to kernel matrix

%------------------------------------------------------
% some algorithms paramaters
%------------------------------------------------------
options.firstbasevariable='first'; % tie breaking method for choosing the base
% variable in the reduced gradient method
options.nbitermax=50;             % maximal number of iteration
options.seuil=0;                   % forcing to zero weights lower than this
options.seuilitermax=10;           % value, for iterations lower than this one

options.miniter=0;                 % minimal number of iterations
options.verbosesvm=0;              % verbosity of inner svm algorithm
options.efficientkernel=1;         % use efficient storage of kernels

err=0;
nsplit= numel(Splits);
dim=size(Kernel);



for i=1:nsplit
    TrainKernel=zeros(size(Splits{i}.indTrain,1),size(Splits{i}.indTrain,1),3);
    TestKernel=zeros(size(Splits{i}.indTest,1),size(Splits{i}.indTrain,1),3);
   
    Kernel(1).stmp(Splits{i}.indTrain,Splits{i}.indTrain);
    TrainKernel(:,:,1)=Kernel(1).stmp(Splits{i}.indTrain,Splits{i}.indTrain);
    TestKernel(:,:,1)=Kernel(1).stmp(Splits{i}.indTest,Splits{i}.indTrain);
    
    Kernel(1).dstmp(Splits{i}.indTrain,Splits{i}.indTrain);
    TrainKernel(:,:,2)=Kernel(1).dstmp(Splits{i}.indTrain,Splits{i}.indTrain);
    TestKernel(:,:,2)=Kernel(1).dstmp(Splits{i}.indTest,Splits{i}.indTrain);
    
    Kernel(1).dstmp(Splits{i}.indTrain,Splits{i}.indTrain);
    TrainKernel(:,:,3)=Kernel(1).ddstmp(Splits{i}.indTrain,Splits{i}.indTrain);
    TestKernel(:,:,3)=Kernel(1).ddstmp(Splits{i}.indTest,Splits{i}.indTrain);
    
    
    tic
    [beta,w,b,posw,story(i),obj(i)] = mklsvm(TrainKernel,T(Splits{i}.indTrain),C,options,verbose);
    
    
    
    ypred=0;
    K1 = TestKernel(:,posw,1);
    K2 = TestKernel(:,posw,2);
    K3 = TestKernel(:,posw,3);
    
    MK = K1*beta(1)+K2*beta(2)+K3*beta(3);
    
    ypred=MK*w+b;
    
    ytest=T(Splits{i}.indTest);
    temp=mean(sign(ypred)==ytest);
    err(i) = temp;

end;

end

