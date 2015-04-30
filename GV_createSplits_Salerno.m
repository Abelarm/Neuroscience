function [splits] = GV_createSplits_Salerno(Labels,params)
% Multivariate decoding of single button press
% project on hemodynamic modeling for Multi-Voxel Pattern Analysis (MVPA)
% Giancarlo Valente, Maastricht, 2014
% giancarlo.valente@maastrichtuniversity.nl


C1                  = [];
for ind             = 1:length(params.classes.C1)
    C1              = [C1 ; find(Labels(:,1) == params.classes.C1(ind))];
end
C2                  = [];
for ind             = 1:length(params.classes.C2)
    C2              = [C2 ; find(Labels(:,1) == params.classes.C2(ind))];
end

Ctot                = union(C1,C2);


switch lower(params.cv.scheme)
    
    case 'lro' % leave run out
        num_runs    = length(unique(params.data.Intervals(:,4)));
        splits      = cell(1,num_runs);
        for ind     = 1:num_runs
            ind_test = find(Labels(:,2)==ind);
            splits{ind}.indTest = ind_test;
            splits{ind}.indTrain = setdiff(Ctot,splits{ind}.indTest);
        end
    case 'kfold' % stratified cross-validation (per class)
        NumFolds    = params.cv.numfolds;
        NumRep      = params.cv.niter;
        ClassSamp   = min(numel(C1),numel(C2)); 
        SampPerFold = floor(ClassSamp/NumFolds);
        splits      = cell(1,NumFolds*NumRep);
        for indRep  = 1:NumRep
            C1p         = C1(randperm(numel(C1)));
            C2p         = C2(randperm(numel(C2)));
            for indFold = 1:NumFolds
                
                if indFold < NumFolds
                    ind_test = 1+SampPerFold*(indFold-1):SampPerFold*indFold;
                else
                    %                     ind_test  = 1:SampPerFold*(indFold-1):numel(C1p);
                    ind_test  = 1+SampPerFold*(indFold-1):ClassSamp;
                end
                C1Test   = C1p(ind_test);
                C1Train  = setdiff(C1p,C1Test);
                
                
                if indFold < NumFolds
                    ind_test = 1+SampPerFold*(indFold-1):SampPerFold*indFold;
                else
%                     ind_test  = 1:SampPerFold*(indFold-1):numel(C2p);
                    ind_test  = 1+SampPerFold*(indFold-1):ClassSamp;
                end
                C2Test   = C2p(ind_test);
                C2Train  = setdiff(C2p,C2Test);
                
                
                splits{NumFolds*(indRep-1) + indFold}.indTrain  = [C1Train(:);C2Train(:)];
                splits{NumFolds*(indRep-1) + indFold}.indTest   = [C1Test(:);C2Test(:)];
            end
        end
    case 'loo'
        splits = cell(1,numel(Ctot));
        for ind = 1:numel(Ctot)
            splits{ind}.indTest    = Ctot(ind);
            splits{ind}.indTrain   = setdiff(Ctot,ind);
        end
        
        
end
       
