function PostAnalysisPermutationTestForLRHemiGLM(tableResults)
    SummaryZResults = [];
    SummaryPvResults = [];
    SummaryPermValueResults = [];
    SummaryRealValueResults = [];
    
    for i = 1:size(tableResults,1)
        cluster1File = [tableResults.GLMVectorLocation{i}, '\AVSD*_Depth2_*.csv'];
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
        SummaryZResults(i) = zscoreForTest;
        SummaryPvResults(i) = pValueForTest;
        SummaryPermValueResults(i,1:permCount) = permutationMeanbetweenHemitree;
        SummaryRealValueResults(i) = realMeanbetweenHemitree;
    end

end