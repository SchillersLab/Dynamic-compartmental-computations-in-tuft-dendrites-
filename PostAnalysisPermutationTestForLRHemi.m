function PostAnalysisPermutationTestForLRHemi(tableResults)
    SummaryZResults = [];
    SummaryPvResults = [];
    SummaryPermValueResults = [];
    SummaryRealValueResults = [];
    
    for i = 1:size(tableResults,1)
        ci = 1;
%         for ci = 1:4
            cluster1File = [tableResults.RunLocation{i}, sprintf('\\All\\ByH\\BetweenAndWithinSubTrees\\AVSD*_Depth2_*.csv')];
            
%             cluster1File = [tableResults.RunLocation{i}, sprintf('\\cluster%d\\ByH\\BetweenAndWithinSubTrees\\AVSD*_Depth2_*.csv', ci)];
            fileList = dir(cluster1File);

            compareValue = [];    
            indexLabel = [];

            for j = 1:size(fileList, 1)
                tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
                if contains(fileList(j).name, 'BetweenMain')
                    indexLabel(end+1:end+length(tableR.Var2)) = 1;
                else
                    indexLabel(end+1:end+length(tableR.Var2)) = 0;
                end

                compareValue(end+1:end+length(tableR.Var2)) = tableR.Var2;
            end

            nanValuesIndex = find(isnan(compareValue));
            
            compareValue(nanValuesIndex) = [];
            indexLabel(nanValuesIndex) = [];
            
            realMeanbetweenHemitree = mean(compareValue(indexLabel==0)) - mean(compareValue(indexLabel==1));

            permCount = 1000;
            permutationMeanbetweenHemitree = nan(1, permCount);
            for j = 1:permCount
                indexP1 = randperm(length(compareValue));
                permutationMeanbetweenHemitree(j) = mean(compareValue(indexLabel(indexP1)==0)) - mean(compareValue(indexLabel(indexP1)==1));
            end

            pValueForTest = sum(permutationMeanbetweenHemitree > realMeanbetweenHemitree)./permCount;

            zscoreForTest = (realMeanbetweenHemitree - mean(permutationMeanbetweenHemitree)) ./ std(permutationMeanbetweenHemitree);
            SummaryZResults(i, ci) = zscoreForTest;
            SummaryPvResults(i,ci) = pValueForTest;
            SummaryPermValueResults(i,ci,1:permCount) = permutationMeanbetweenHemitree;
            SummaryRealValueResults(i,ci) = realMeanbetweenHemitree;
%         end
    end
end