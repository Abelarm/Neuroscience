% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl


clear variables;
close all;

%FullPath= '/Users/Luigi/Desktop/Dropbox/Personali/UNIVERSITA/IC/ProgettoIC/progetto_neuroscienze/Data/';
FullPath='D:\Desktop\Dropbox\Progetti IC 2014\progetto neuroscienze\Data';
Subject = 'JE'; % it can be 'AZ', 'GV' or 'JE';

% add Neuroelf to the path (change the path accordingly!
%addpath(genpath('/Users/Luigi/Desktop/Dropbox/Personali/UNIVERSITA/IC/ProgettoIC/progetto_neuroscienze/Toolboxes/NeuroElf_v10_4562'))
addpath(genpath('D:\Desktop\Dropbox\ProgettoIC\progetto_neuroscienze\Toolboxes\NeuroElf_v10_4562'));
% adding Spider to the path (change the path accordingly!)
%addpath(genpath('/Users/Luigi/Desktop/Dropbox/Personali/UNIVERSITA/IC/ProgettoIC/progetto_neuroscienze/Toolboxes/spider'));
addpath(genpath('D:\Desktop\Dropbox\ProgettoIC\progetto_neuroscienze\Toolboxes\spider'));

%try
    % some initialization: getting the names of files, masks and the
    % parameters.
    %NameFiles                       =  GV_getNameFiles_Salerno(Subject);
    %Params                          =  GV_initializeParams_Salerno;
    
    
    % read data and protocols
    %[Data,Params]                   = GV_importfMRIdata_Salerno(NameFiles,Params);
    %[DM,Labels,Params]              = GV_importProtocol_Salerno(NameFiles,Params);
    % extract Trials
    %[Trials]                        = GV_getTrials_Salerno(Data,Labels,Params);
%catch
    % if for some reason the BV reader does not work or the data are not
    % accessible, all the data are here
    load([FullPath '\' Subject '\' Subject '_data.mat']);
%end
% define splitting scheme
[Trials,Labels]                 = GR_CleanTrials(Trials,Labels);
[Splits]                        = GV_createSplits_Salerno(Labels,Params);
% feature extraction, you will work here!
[X]                             = GV_extractFeatures_Salerno(Trials,Splits,Params);

[TrialsMed]                     = GR_medTrials(Trials);
[Data8Kernel]                   = GR_genPointDataKernel(TrialsMed);

[HRFBasis]                      = GR_extractHRFBasis_Salerno();
[Data3Kernel]                   = GR_extractFeatures3KernelData_Salerno(TrialsMed,HRFBasis);

[KernelPoints]                  = GR_gen8LinearKernel(Data8Kernel);
[KernelBasis]                   = GR_gen3LinearKernel(Data3Kernel);
[T]                             = GR_getT(Trials);

NITER=30;
%tic
%ERR3=GR_nedstedCrossValidation(KernelBasis,T,Splits);
%toc
% [ ERR8 ]                         = GR_classifyCrossVal8Kernel( KernelPoints,Splits,T,NITER);
 Parame(2)=3;
 Parame(1)=50;
[ ERR3 ]                         = GR_classifyCrossVal3Kernel( KernelBasis,Splits,T,Parame,NITER);
%[ ERR11 ]                        = GR_classifyCrossVal11Kernel( KernelPoints,KernelBasis,Splits,T,Parame,NITER);
  
%[ ERR1 ]                          =GR_SimpleMKL3(KernelBasis,T,Splits);
%[ ERR2 ]                          =GR_SimpleMKL8(KernelPoints,T,Splits);

%medERR1                           =GR_CalcolaMedia(ERR1,size(Trials,2))
%medERR2                           =GR_CalcolaMedia(ERR2,size(Trials,2))
medERR3                           =GR_CalcolaMedia(ERR3,size(Trials,2))


% classify
[Results]                         = GV_classifyCrossVal_Salerno(X,Labels,Splits);
%ERR(1:50)                         =  Results.mean_err;

%plot(ERR,'g');
%hold all
%plot(ERR1,'b');
%plot(ERR2,'r');
%plot(ERR3,'k');

disp(Results)


MaskSize                        = [87 60 69];
%[NameSave]                      = GV_writemapMVPA_Salerno(Params,NameFiles,Results,MaskSize);
%save(NameSave)