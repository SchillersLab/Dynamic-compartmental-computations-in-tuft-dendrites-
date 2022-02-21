function PermutationTestForLRHemi(roiActivityDistanceMatrixByH,selectedROI,selectedROISplitDepth1, outputpath)
    compareValue = [];
    indexLabel = [];
    for i = 1:length(selectedROI)
        for j = (i+1):length(selectedROI)
            if isnan(roiActivityDistanceMatrixByH(i,j))
                continue;
            end
            
            compareValue(end+1) = roiActivityDistanceMatrixByH(i,j);
            indexLabel(end+1) = selectedROISplitDepth1(i) ~= selectedROISplitDepth1(j);
        end
    end
    
%     realMeanbetweenHemitree = mean(compareValue(indexLabel==1));
    realMeanbetweenHemitree = mean(compareValue(indexLabel==0)) - mean(compareValue(indexLabel==1));
    
    permCount = 1000;
    permutationMeanbetweenHemitree = nan(1, permCount);
    for i = 1:permCount
        indexP1 = randperm(length(compareValue));
%         permutationMeanbetweenHemitree(i) = mean(compareValue(indexLabel(indexP1)==1));
        permutationMeanbetweenHemitree(i) = mean(compareValue(indexLabel(indexP1)==0)) - mean(compareValue(indexLabel(indexP1)==1));
    end
    
%     pValueForTest = sum(permutationMeanbetweenHemitree < realMeanbetweenHemitree)./permCount;
    pValueForTest = sum(permutationMeanbetweenHemitree > realMeanbetweenHemitree)./permCount;
    
    zscoreForTest = (realMeanbetweenHemitree - mean(permutationMeanbetweenHemitree)) ./ std(permutationMeanbetweenHemitree);
    
%     save([outputpath, '\permutationTest_glmVector.mat'], 'pValueForTest', 'zscoreForTest', 'permutationMeanbetweenHemitree', 'realMeanbetweenHemitree');
        save([outputpath, '\permutationTestWithinVsBetween_Vector.mat'], 'pValueForTest', 'zscoreForTest', 'permutationMeanbetweenHemitree', 'realMeanbetweenHemitree');
end