function [results] = GV_classifyCrossVal_Salerno(X,Labels,Splits)
% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl

Labelsgroup         = Labels(:,1);
Labelsgroup(Labelsgroup==4) = 0;


UniqueLabels        =  setdiff(unique(Labelsgroup),0);

err                 = zeros(numel(Splits),1);
err1                = zeros(numel(Splits),1);
err2                = zeros(numel(Splits),1);
test_size           = zeros(numel(Splits),1);
test_size1          = zeros(numel(Splits),1);
test_size2          = zeros(numel(Splits),1);
temp_model_weights  = zeros(size(X,2),length(Splits));

for ind= 1: numel(Splits)
%     fprintf('.');
    if size(X,3)==1
        Data_train  = X(Splits{ind}.indTrain,:);
        Data_test   = X(Splits{ind}.indTest,:);
    else
        Data_train  = X(Splits{ind}.indTrain,:,ind);
        Data_test   = X(Splits{ind}.indTest,:,ind);
    end
    Labels_train    = Labelsgroup(Splits{ind}.indTrain);
    Labels_test     = Labelsgroup(Splits{ind}.indTest);
    
    
    %esegue 50 volte SVM sui diversi split
    [res] = GV_traintestSVM_Salerno(Data_train,Data_test,Labels_train,Labels_test);
   
    
    err(ind) = res.err;
    err1(ind)= res.errs1;
    err2(ind)= res.errs2;
    test_size(ind) = size(Data_test,1);
    test_size1(ind) = sum(Labels_test == UniqueLabels(1));
    test_size2(ind) = sum(Labels_test == UniqueLabels(2));
    temp_model_weights(:,ind) = res.model_weights(:);
    
    
end
%calcola gli errori totali per ogni iterazione
results.err = err;
results.err1= err1;
results.err2= err2;
results.tot_err = sum(err);
results.tot_err1= sum(err1);
results.tot_err2= sum(err2);
results.mean_err= sum(err)/sum(test_size);
results.mean_err1 = sum(err1)/sum(test_size1);
results.mean_err2 = sum(err2)/sum(test_size2);
results.temp_model_weights = temp_model_weights;
% fprintf('\n');
