function [ splitsN ] = GR_split88_72(  indexSplit )
%   Scrivo l'array ricevuto in input, in ogni riga della matrice finale
%   e poi effettuo lo split
    
    sizeSplit=size(indexSplit,1);
    if(sizeSplit==88)
        tmp=1;
        splitM=zeros(11,88);
        splitsN= cell(1,11);
        for i=1:11
            splitM(i,:)=indexSplit;
        end
        splitsN{1}.indTrain = splitM(1,1:80);
        splitsN{1}.indTest = splitM(1,81:88);
    %   Adesso scambio, in modo che i primi 80 di ogni riga sono per training e gli ultimi 8    
        for i=2:11      
              splitM(i,[tmp:(tmp+7) 81:88])=splitM(i,[81:88 tmp:(tmp+7)]);
            splitsN{i}.indTrain = splitM(i,1:80);
            splitsN{i}.indTest = splitM(i,81:88);
            tmp=tmp+8;
        end 
    end
    if(sizeSplit==72)
        tmp=1;
        splitM=zeros(9,72);
        splitsN= cell(1,9);
        for i=1:9
            splitM(i,:)=indexSplit;
        end
        splitsN{1}.indTrain = splitM(1,1:64);
        splitsN{1}.indTest = splitM(1,65:72);
    %   Adesso scambio, in modo che i primi 64 di ogni riga sono per training e gli ultimi 8    
        for i=2:9      
              splitM(i,[tmp:(tmp+7) 65:72])=splitM(i,[65:72 tmp:(tmp+7)]);
            splitsN{i}.indTrain = splitM(i,1:64);
            splitsN{i}.indTest = splitM(i,65:72);
            tmp=tmp+8;
        end 
    end
            
end

