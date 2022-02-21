function plotMantelResults(clusterCount, eventTable, outputpathCurr, clusterType, otherPath)
    clusterP = cell(clusterCount + 1, 2);
    for i = 0:clusterCount
        curP = eventTable.MantelPval(eventTable.clusterByH == i | i == 0);
        clusterP(i+1, 1) = {sum(curP <= 0.05) ./ length(curP)};
        clusterP(i+1, 2) = {sprintf('cluster_%d', i)};
    end
    
    clusterP(1,2) = {replace(clusterP{1,2}, 'cluster_0', 'All')};
    
    fig = figure;
    hold on;
    title('Mantel Test Results');
    xlabel('#cluster (By High)');
    ylabel('Precentage of Events Significant correlated to structure');
    bar(0:clusterCount, cell2mat(clusterP(:,1)));
    xlim([-1,clusterCount+1]);
    xticks(0:clusterCount);
    xticklabels(clusterP(:,2));
    
    mysave(fig, fullfile(outputpathCurr, 'mentalResultsSummary'));
    
    writetable(cell2table(clusterP), fullfile(outputpathCurr, 'mentalResultsSummary.csv'));
    
    mysave(fig, fullfile(otherPath, 'mentalResultsSummary'));
    
    writetable(cell2table(clusterP), fullfile(otherPath, 'mentalResultsSummary.csv'));
end