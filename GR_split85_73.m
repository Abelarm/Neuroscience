function [ splitsN ] = split88_72(  indexSplit )
%   Scrivo l'array ricevuto in input, in ogni riga della matrice finale
%   e poi effettuo lo split
    
    sizeSplit=size(indexSplit,1);
    if(sizeSplit==85)
        tmp=1;
        splitM=zeros(17,85);
        splitsN= cell(1,17);
        for i=1:17
            splitM(i,:)=indexSplit;
        end
        splitsN{1}.indTrain = splitM(1,1:80);
        splitsN{1}.indTest = splitM(1,81:85);
    %   Adesso scambio, in modo che i primi 80 di ogni riga sono per training e gli ultimi 5 test 
        for i=2:17      
              splitM(i,[tmp:(tmp+4) 81:85])=splitM(i,[81:85 tmp:(tmp+4)]);
            splitsN{i}.indTrain = splitM(i,1:80);
            splitsN{i}.indTest = splitM(i,81:85);
            tmp=tmp+5;
        end 
    end
    if(sizeSplit==73)
        tmp=1;
        splitM=zeros(73,73);
        splitsN= cell(1,73);
        for i=1:73
            splitM(i,:)=indexSplit;
        end
        splitsN{1}.indTrain = splitM(1,1:72);
        splitsN{1}.indTest = splitM(1,73);
    %   Adesso scambio, in modo che i primi 72 di ogni riga sono per training e l'ultimo 1 test    
        for i=2:73      
              splitM(i,[tmp 73])=splitM(i,[73 tmp]);
            splitsN{i}.indTrain = splitM(i,1:72);
            splitsN{i}.indTest = splitM(i,73);
            tmp=tmp+1;
        end 
    end         
end

