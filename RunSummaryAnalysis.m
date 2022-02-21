function RunSummaryAnalysis
    mantelRunnerLocation = '\\jackie-analysis\e\Shay\RunnersLocationSummary.xlsx';
    sheetName = 'RunOnlyTuft';
    outputfolder = '\\jackie-analysis\e\Shay\MantelPostSummary\Big_19_10\';
    mkdir(outputfolder);
    classification = 1;
    
    tableResults = readtable(mantelRunnerLocation,'Sheet',sheetName);
    summaryTMantel = zeros(4, sum(tableResults.Classification == classification & tableResults.includeMainR2 == 1));
    summaryTDist = zeros(4, sum(tableResults.Classification == classification & tableResults.includeMainR2 == 1));
    summaryMantel = zeros(4, sum(tableResults.Classification == classification & tableResults.includeMainR2 == 1));
    
    indexLocation = find(tableResults.Classification == classification & tableResults.includeMainR2 == 1);
    
%     summaryEventsHistogram = zeros(4, sum(tableResults.Classification == classification));
    for i = 1:size(indexLocation,1)
        
%         eventsClusterSummary = readtable(fullfile(tableResults.RunLocation{indexLocation(i)}, 'HistogramEventsClusterSummary.csv'));
%         
%         col1 = contains(eventsClusterSummary.Cluster, '1');
%         col2 = contains(eventsClusterSummary.Cluster, '2');
%         col3 = contains(eventsClusterSummary.Cluster, '3');
%         col4 = contains(eventsClusterSummary.Cluster, '4');
%         
%         summaryEventsHistogram(1, i) = eventsClusterSummary.ByHVAlue(col1);
%         summaryEventsHistogram(2, i) = eventsClusterSummary.ByHVAlue(col2);
%         summaryEventsHistogram(3, i) = eventsClusterSummary.ByHVAlue(col3);
%         summaryEventsHistogram(4, i) = eventsClusterSummary.ByHVAlue(col4); 
%         
        mantelbycluster = readtable(fullfile(tableResults.RunLocation{indexLocation(i)}, '\Mantel\Reg\PerClusterMantelTable.csv'));
        
        colAll = contains(mantelbycluster.Cluster, 'All');
        col1 = contains(mantelbycluster.Cluster, '1');
        col2 = contains(mantelbycluster.Cluster, '2');
        col3 = contains(mantelbycluster.Cluster, '3');
        col4 = contains(mantelbycluster.Cluster, '4');
        
        summaryMantel(1, i) = mantelbycluster.MantelValue(col1);
        summaryMantel(2, i) = mantelbycluster.MantelValue(col2);
        summaryMantel(3, i) = mantelbycluster.MantelValue(col3);
        summaryMantel(4, i) = mantelbycluster.MantelValue(col4);
        summaryMantel(5, i) = mantelbycluster.MantelValue(colAll);
        
        mT1 = load(fullfile(tableResults.RunLocation{indexLocation(i)}, '\Mantel\Threshold\ShortestPath_cluster1.mat'));
        mT2 = load(fullfile(tableResults.RunLocation{indexLocation(i)}, '\Mantel\Threshold\ShortestPath_cluster2.mat'));
        mT3 = load(fullfile(tableResults.RunLocation{indexLocation(i)}, '\Mantel\Threshold\ShortestPath_cluster3.mat'));
        mT4 = load(fullfile(tableResults.RunLocation{indexLocation(i)}, '\Mantel\Threshold\ShortestPath_cluster4.mat'));
        
        summaryTMantel(1, i) = mT1.maxMantel;
        summaryTMantel(2, i) = mT2.maxMantel;
        summaryTMantel(3, i) = mT3.maxMantel;
        summaryTMantel(4, i) = mT4.maxMantel;
        
        
        if isfield(mT1, 'tJump')
            test = (length(mT4.outM)-1) * mT1.tJump;
        else
            test = (length(mT4.outM)-1) * 50;
        end
        
        if isnan(mT1.maxThreshold)
            summaryTDist(1, i) = 0;
        else
            summaryTDist(1, i) = mT1.maxThreshold ./ test;
        end
        
        if isnan(mT2.maxThreshold)
            summaryTDist(2, i) = 0;
        else
            summaryTDist(2, i) = mT2.maxThreshold ./ test;
        end
        
        if isnan(mT4.maxThreshold)
            summaryTDist(4, i) = 0;
        else
            summaryTDist(4, i) = mT4.maxThreshold ./ test;
        end
        
        
        if isnan(mT3.maxThreshold)
            summaryTDist(3, i) = 0;
        else
            summaryTDist(3, i) = mT3.maxThreshold ./ test;
        end
       
    end    

    namePlot = {'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4'};
   
%% Section - statistic and Ploting Threshold dist
    thresholdF = figure;
    hold on;
    
    boxchart(summaryTDist');
    
    xticklabels(namePlot);
    xtickangle(90);
    ylabel('Selected Threshold (microns) Per Cluster');
    
    meanTM = mean(summaryTDist, 2);
    plot(meanTM,'-*');
    
    textAnova = '';
    [p,~,statsTDist] = anova1(summaryTDist',namePlot, 'off');
    if p < 0.05
        [c,~,~,groupnames] = multcompare(statsTDist, 'Display', 'off');
        for j = 1:size(c, 1)
            textAnova = strcat(textAnova, sprintf('%s vs %s - pValue : %.4f, mean Diff : %.4f, CI_L : %.4f, CI_H : %.4f \\n ', ...
                groupnames{c(j, 1)}, groupnames{c(j, 2)}, c(j, 6), c(j, 4), c(j, 3), c(j, 5)));
        end
    end
    
    fid=fopen(fullfile(outputfolder, 'statistic_MantelDistThreshold.txt'),'w');
    fprintf(fid, textAnova);
    fclose(fid);
    
    mysave(thresholdF, fullfile(outputfolder, 'MantelDistThresholdSummaryPlot'));   
    
    textStat = '';
    for j = 1:size(summaryTDist, 1)
        stdTM(j) = std(summaryTDist(j, :));
        textStat = strcat(textStat, sprintf('cluster %d, mean = %f, std = %f',j, meanTM(j), stdTM(j)));
    end
    
    fid=fopen(fullfile(outputfolder, 'statisticMeanSTD_MantelDistThreshold.txt'),'w');
    fprintf(fid, textStat);
    fclose(fid);
    
    
%% Section - statistic and Ploting Threshold Mantel
    mantelF = figure;
    hold on;
    
    boxchart(summaryTMantel');
    xticklabels(namePlot);
    xtickangle(90);
    ylabel('Selected Threshold - Mantel Value Per Cluster');

    meanTM = mean(summaryTMantel, 2);
    plot(meanTM,'-*');
    
    textAnova = '';
    [p,~,statsTM] = anova1(summaryTMantel', namePlot, 'off');
    if p < 0.05
        [c,~,~,groupnames] = multcompare(statsTM, 'Display', 'off');
        for j = 1:size(c, 1)
            textAnova = strcat(textAnova, sprintf('%s vs %s - pValue : %.4f, mean Diff : %.4f, CI_L : %.4f, CI_H : %.4f \\n ', ...
                groupnames{c(j, 1)}, groupnames{c(j, 2)}, c(j, 6), c(j, 4), c(j, 3), c(j, 5)));
        end
    end
    
    fid=fopen(fullfile(outputfolder, 'statistic_MantelVThreshold.txt'),'w');
    fprintf(fid, textAnova);
    fclose(fid);
    
    mysave(mantelF, fullfile(outputfolder, 'MantelThresholdSummaryPlot'));
    
    textStat = '';
    for j = 1:size(summaryTMantel, 1)
        stdTM(j) = std(summaryTMantel(j, :));
        textStat = strcat(textStat, sprintf('cluster %d, mean = %f, std = %f',j, meanTM(j), stdTM(j)));
    end
    
    fid=fopen(fullfile(outputfolder, 'statisticMeanSTD_MantelVThreshold.txt'),'w');
    fprintf(fid, textStat);
    fclose(fid);
   
    
%% Section - statistic and Ploting Mantel Per Cluster
    mantelPerClusterF = figure;
    hold on;
    
    boxchart(summaryMantel');
    xticklabels(namePlot);
    xtickangle(90);
    ylabel('Mantel Value Per Cluster');

    meanTM = mean(summaryMantel, 2);
    plot(meanTM,'-*');
    
    textAnova = '';
    [p,~,statsM] = anova1(summaryMantel', namePlot, 'off');
    if p < 0.05
        [c,~,~,groupnames] = multcompare(statsM, 'Display', 'off');
        for j = 1:size(c, 1)
            textAnova = strcat(textAnova, sprintf('%s vs %s - pValue : %.4f, mean Diff : %.4f, CI_L : %.4f, CI_H : %.4f \\n ', ...
                groupnames{c(j, 1)}, groupnames{c(j, 2)}, c(j, 6), c(j, 4), c(j, 3), c(j, 5)));
        end
    end
    
    fid=fopen(fullfile(outputfolder, 'statistic_MantelVPerCluster.txt'),'w');
    fprintf(fid, textAnova);
    fclose(fid);
    
    mysave(mantelPerClusterF, fullfile(outputfolder, 'MantelPerClusterSummaryPlot'));   
    
    
    textStat = '';
    for j = 1:size(summaryMantel, 1)
        stdTM(j) = std(summaryMantel(j, :));
        textStat = strcat(textStat, sprintf('cluster %d, mean = %f, std = %f',j, meanTM(j), stdTM(j)));
    end
    
    fid=fopen(fullfile(outputfolder, 'statisticMeanSTD_MantelPerClusterSummary.txt'),'w');
    fprintf(fid, textStat);
    fclose(fid);
   
end