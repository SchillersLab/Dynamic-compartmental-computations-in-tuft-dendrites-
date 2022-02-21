function PostAnalysisPermutationTestForR2(tableResults)
    SummaryZResults = [];
    SummaryPvResults = [];
    SummaryPermValueResults = [];
    SummaryRealValueResults = [];
    
    for i = 1:size(tableResults,1)
        for ci = 1:4
            cluster1File = [tableResults.RunLocation{i}, sprintf('\\cluster%d\\ByH\\BetweenAndWithinSubTrees\\AVSD*_Depth2_*.csv', ci)];
            fileList = dir(cluster1File);
            X_value = [];
            Y_value = [];
            for j = 1:size(fileList, 1)
                tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
                X_value(end+1:end+length(tableR.Var2)) = tableR.Var1;
                Y_value(end+1:end+length(tableR.Var2)) = tableR.Var2;
            end
            
            mdAll = fitglm(X_value, Y_value);
            trueR2_value = mdAll.Rsquared.Ordinary;

            [pValueForTest, permutationR2]  = foldsCalculator(X_value,Y_value, trueR2_value);
            
            zscoreForTest = (trueR2_value - mean(permutationR2)) ./ std(permutationR2);
            SummaryZResults(i, ci) = zscoreForTest;
            SummaryPvResults(i,ci) = pValueForTest;
            SummaryPermValueResults(i,ci,1:length(permutationR2)) = permutationR2;
            SummaryRealValueResults(i,ci) = trueR2_value;
        end
    end
end