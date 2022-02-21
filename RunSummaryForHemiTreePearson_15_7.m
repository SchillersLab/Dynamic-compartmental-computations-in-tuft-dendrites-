function RunSummaryForHemiTreePearson_15_7()
    mantelRunnerLocation = '\\jackie-analysis\e\Shay\RunnersLocationSummary.xlsx';
    sheetName = 'RunOnlyTuft';
    outputfolder(1) = {'\\jackie-analysis\e\Shay\StatisticSummary\ByH\Big_2002_Final1\'};
    mkdir(outputfolder{1});
    classification(1) = 1;
    
    outputfolder(2) = {'\\jackie-analysis\e\Shay\StatisticSummary\ByH\Small_2002_Final1\'};
    mkdir(outputfolder{2});
    classification(2) = 2;
    
    colorsPerType = zeros([2,3]);
    colorsPerType(2, :) = [0,0,1];
    
    pcaMeanAll = {};
    pcaChance = {};
    pcaSTDAll = {};
    matAll = [];
    matAll2 = [];
    slopAll = {};
    tableResults = readtable(mantelRunnerLocation,'Sheet',sheetName);
    indexLocation1 = find(tableResults.Classification == classification(1) & tableResults.R2HemiTreeInclude == 1);
    indexLocation2 = find(tableResults.Classification == classification(2) & tableResults.R2HemiTreeInclude == 1);
    
    indexLocationR2_1 = find(tableResults.Classification == classification(1) & tableResults.includeMainR2 == 1);
    indexLocationR2_2 = find(tableResults.Classification == classification(2) & tableResults.includeMainR2 == 1);
    
%     h_all = [];
%     for i = 1:size(tableResults, 1)
%         events = readtable(fullfile(tableResults.RunLocation{i}, 'eventsCaSummary.csv'));
%         h_all(end+1:end+size(events,1), 1:2) = [events.clusterByH,events.H];
%     end
%     figure;hold on; 
%     histfit(h_all(:, 2),53, 'lognormal')
%     
    for k = 1:2
        indexLocation = find(tableResults.Classification == classification(k) & tableResults.R2HemiTreeInclude == 1);
        indexLocationMain = find(tableResults.Classification == classification(k) & tableResults.includeMainR2 == 1);
        indexLocationCompare = find(tableResults.Classification == classification(k) & tableResults.includeMainPersonCompare2 == 1);
       
        
        matRTemp2 = BetweenVsWithinHemiDepth(indexLocationCompare, tableResults, outputfolder{k}, 0.05);
        matRTemp2(:, 3) = k;
        
        matAll2(end+1:end + size(matRTemp2, 1),1:4) = matRTemp2;
        
        
        [matRTemp, bslop] = BetweenVsWithinR2(indexLocationMain,indexLocation, tableResults, outputfolder{k}, 0.05);
        matRTemp(:, 3) = k;
        slopAll(k) = {bslop};

        matAll(end+1:end + size(matRTemp, 1),1:4) = matRTemp;
        
        statPlot([tableResults.R2Cluster1_ALL(indexLocationMain),...
            tableResults.R2Cluster2_ALL(indexLocationMain),...
            tableResults.R2Cluster3_ALL(indexLocationMain),...
            tableResults.R2Cluster4_ALL(indexLocationMain)], 'R2', 'R2Main', outputfolder{k}, colorsPerType(k, :))
          
        statPlot([tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation),...
            tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation),...
            tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation),...
            tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation)], 'R2', 'R2maxHemiTree', outputfolder{k}, colorsPerType(k, :))

        statPlot([tableResults.R2Cluster1_HemiTree(indexLocation),...
            tableResults.R2Cluster2_HemiTree(indexLocation),...
            tableResults.R2Cluster3_HemiTree(indexLocation),...
            tableResults.R2Cluster4_HemiTree(indexLocation)], 'R2', 'R2AllHemiTree', outputfolder{k}, colorsPerType(k, :))

        pcaMean = [];
        pcaSTD = [];
        for i = 1:size(indexLocationMain,1)
            spR = split(tableResults.PcaAccuracyMeanSTD_cluster1{indexLocationMain(i)}, ',');
            pcaMean(i,1) = str2double(spR{1});        
            pcaSTD(i,1) = str2double(spR{2});
            spR = split(tableResults.PcaAccuracyMeanSTD_cluster2{indexLocationMain(i)}, ',');
            pcaMean(i,2) = str2double(spR{1});
            pcaSTD(i,2) = str2double(spR{2});
            spR = split(tableResults.PcaAccuracyMeanSTD_cluster3{indexLocationMain(i)}, ',');
            pcaMean(i,3) = str2double(spR{1});
            pcaSTD(i,3) = str2double(spR{2});
            spR = split(tableResults.PcaAccuracyMeanSTD_cluster4{indexLocationMain(i)}, ',');
            pcaMean(i,4) = str2double(spR{1});
            pcaSTD(i,4) = str2double(spR{2});
        end

        statPlot([pcaMean(:,1),...
            pcaMean(:,2),...
            pcaMean(:,3),...
            pcaMean(:,4)], 'PCA Accuracy', 'PcaAccuracy', outputfolder{k}, colorsPerType(k, :))

        statPlot([pcaSTD(:,1),...
            pcaSTD(:,2),...
            pcaSTD(:,3),...
            pcaSTD(:,4)], 'PCA STD Accuracy', 'PcaAccuracySTD', outputfolder{k}, colorsPerType(k, :))

        
        for j = 1:4
            f = figure;hold on;
            shadedErrorBar([], pcaMean(:,j), pcaSTD(:,j), 'lineProps', '-k');
            plot(1:size(pcaMean,1), tableResults.PcaChanceLevel(indexLocationMain), '--k', 'LineWidth', 2);
            title(sprintf('Cluster %d, Pca Mean + STD + Chance ', j));
            xlabel('#Neuron');
            ylabel('PCA Accuracy');
            mysave(f, fullfile(outputfolder{k}, sprintf('PCASummary_cluster%d', j)));
        end
        
        f = figure;hold on;
        title('PCA accuracy Mean and STD')
        errorbar(mean(pcaMean), mean(pcaSTD));
        plot(1:size(pcaMean, 2), ones(1, size(pcaMean, 2))*mean(tableResults.PcaChanceLevel(indexLocationMain)), '--k', 'LineWidth', 2)
        
        xticks(1:size(pcaMean, 2));
        xticklabels({'cluster1', 'cluster2', 'cluster3', 'cluster4'});
        xlim([0,size(pcaMean, 2)+1]);
        mysave(f, fullfile(outputfolder{k}, 'PCASummary_All'));
        
        pcaMeanAll(k) = {pcaMean};
        pcaSTDAll(k) = {pcaSTD};
        pcaChance(k) = {mean(tableResults.PcaChanceLevel(indexLocationMain))};
        
        text = '';
        text = strcat(text, sprintf('Cluster1 pass chance %f \\n', mean(pcaMean(:, 1) > tableResults.PcaChanceLevel(indexLocationMain))));
        text = strcat(text, sprintf('Cluster2 pass chance %f \\n', mean(pcaMean(:, 2) > tableResults.PcaChanceLevel(indexLocationMain))));
        text = strcat(text, sprintf('Cluster3 pass chance %f \\n', mean(pcaMean(:, 3) > tableResults.PcaChanceLevel(indexLocationMain))));
        text = strcat(text, sprintf('Cluster4 pass chance %f \\n', mean(pcaMean(:, 4) > tableResults.PcaChanceLevel(indexLocationMain))));
        
        fid=fopen(fullfile(outputfolder{k}, 'statisticPCAPassChanceLevel.txt'),'w');
        fprintf(fid, text);
        fclose(fid);

        
        text = '';
        text = strcat(text, sprintf('Cluster1 mean %f \\n', mean(tableResults.AnovaTestBetweenAndWithIn_cluster1(indexLocationMain))));
        text = strcat(text, sprintf('Cluster2 mean %f \\n', mean(tableResults.AnovaTestBetweenAndWithIn_cluster2(indexLocationMain))));
        text = strcat(text, sprintf('Cluster3 mean %f \\n', mean(tableResults.AnovaTestBetweenAndWithIn_cluster3(indexLocationMain))));
        text = strcat(text, sprintf('Cluster4 mean %f \\n', mean(tableResults.AnovaTestBetweenAndWithIn_cluster4(indexLocationMain))));

        fid=fopen(fullfile(outputfolder{k}, 'statisticMeanAnovaSignificantBetweenGroups.txt'),'w');
        fprintf(fid, text);
        fclose(fid);
    end
    
    f = figure;hold on;
    title('Slop')
  
    bar(categorical({'T1_1', 'T1_2', 'T1_3', 'T1_4', 'T2_1', 'T2_2', 'T2_3', 'T2_4'},...
        {'T1_1', 'T2_1', 'T1_2', 'T2_2', 'T1_3', 'T2_3', 'T1_4', 'T2_4'}),[mean(slopAll{1}),mean(slopAll{2})]);
        
    mysave(f, fullfile(outputfolder{1}, 'SlopSummary_All_BigAndSmall'));
    
    [p1, h1] = ranksum(slopAll{1}(:,1),  slopAll{2}(:,1));
    [p2, h2] = ranksum(slopAll{1}(:,2),  slopAll{2}(:,2));
    [p3, h3] = ranksum(slopAll{1}(:,3),  slopAll{2}(:,3));
    [p4, h4] = ranksum(slopAll{1}(:,4),  slopAll{2}(:,4));

    text = '';
    text = strcat(text, sprintf('Cluster1 h1 %f, p %f \\n', h1,p1));
    text = strcat(text, sprintf('Cluster2 h1 %f, p %f \\n', h2,p2));
    text = strcat(text, sprintf('Cluster3 h1 %f, p %f \\n', h3,p3));
    text = strcat(text, sprintf('Cluster4 h1 %f, p %f \\n', h4,p4));

    fid=fopen(fullfile(outputfolder{1}, 'statisticranksumSlopSmallVsBig.txt'),'w');
    fprintf(fid, text);
    fclose(fid);
    
    f = figure;hold on;
    title('R2')
     bar(categorical({'T1_1', 'T1_2', 'T1_3', 'T1_4', 'T2_1', 'T2_2', 'T2_3', 'T2_4'},...
        {'T1_1', 'T2_1', 'T1_2', 'T2_2', 'T1_3', 'T2_3', 'T1_4', 'T2_4'}),[mean(tableResults.R2Cluster1_ALL(indexLocationR2_1)), ...
        mean(tableResults.R2Cluster2_ALL(indexLocationR2_1)),...
        mean(tableResults.R2Cluster3_ALL(indexLocationR2_1)),...
        mean(tableResults.R2Cluster4_ALL(indexLocationR2_1)),...
        mean(tableResults.R2Cluster1_ALL(indexLocationR2_2)), ...
        mean(tableResults.R2Cluster2_ALL(indexLocationR2_2)),...
        mean(tableResults.R2Cluster3_ALL(indexLocationR2_2)),...
        mean(tableResults.R2Cluster4_ALL(indexLocationR2_2))]);
       
    mysave(f, fullfile(outputfolder{1}, 'R2Summary_All_BigAndSmall'));
    
    [p1, h1] = ranksum(tableResults.R2Cluster1_ALL(indexLocationR2_1),  tableResults.R2Cluster1_ALL(indexLocationR2_2));
    [p2, h2] = ranksum(tableResults.R2Cluster2_ALL(indexLocationR2_1),  tableResults.R2Cluster2_ALL(indexLocationR2_2));
    [p3, h3] = ranksum(tableResults.R2Cluster3_ALL(indexLocationR2_1),  tableResults.R2Cluster3_ALL(indexLocationR2_2));
    [p4, h4] = ranksum(tableResults.R2Cluster4_ALL(indexLocationR2_1),  tableResults.R2Cluster4_ALL(indexLocationR2_2));

    text = '';
    text = strcat(text, sprintf('Cluster1 h1 %f, p %f \\n', h1,p1));
    text = strcat(text, sprintf('Cluster2 h1 %f, p %f \\n', h2,p2));
    text = strcat(text, sprintf('Cluster3 h1 %f, p %f \\n', h3,p3));
    text = strcat(text, sprintf('Cluster4 h1 %f, p %f \\n', h4,p4));

    fid=fopen(fullfile(outputfolder{1}, 'statisticranksumR2SmallVsBig.txt'),'w');
    fprintf(fid, text);
    fclose(fid);
    
    f = figure;hold on;
    title('PCA accuracy Mean')
    
    meanAll(1:4) = mean(pcaMeanAll{1});
    meanAll(5:8) = mean(pcaMeanAll{2});
    bar(categorical({'T1_1', 'T1_2', 'T1_3', 'T1_4', 'T2_1', 'T2_2', 'T2_3', 'T2_4'},...
        {'T1_1', 'T2_1', 'T1_2', 'T2_2', 'T1_3', 'T2_3', 'T1_4', 'T2_4'}),meanAll);

    mysave(f, fullfile(outputfolder{1}, 'PCASummary_All_BigAndSmall'));
    
    [p1, h1] = ranksum(pcaMeanAll{1}(:,1),  pcaMeanAll{2}(:,1));
    [p2, h2] = ranksum(pcaMeanAll{1}(:,2),  pcaMeanAll{2}(:,2));
    [p3, h3] = ranksum(pcaMeanAll{1}(:,3),  pcaMeanAll{2}(:,3));
    [p4, h4] = ranksum(pcaMeanAll{1}(:,4), pcaMeanAll{2}(:,4));

    text = '';
    text = strcat(text, sprintf('Cluster1 h1 %f, p %f \\n', h1,p1));
    text = strcat(text, sprintf('Cluster2 h1 %f, p %f \\n', h2,p2));
    text = strcat(text, sprintf('Cluster3 h1 %f, p %f \\n', h3,p3));
    text = strcat(text, sprintf('Cluster4 h1 %f, p %f \\n', h4,p4));

    fid=fopen(fullfile(outputfolder{1}, 'statisticranksumPCASmallVsBig.txt'),'w');
    fprintf(fid, text);
    fclose(fid);
    
    
    f = figure; hold on;

    meanAll(1) = mean(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation1));
    meanAll(2) = mean(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation1));
    meanAll(3) = mean(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation1));   
    meanAll(4) = mean(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation1));  
    meanAll(5) = mean(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation2));
    meanAll(6) = mean(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation2));    
    meanAll(7) = mean(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation2));  
    meanAll(8) = mean(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation2));
     bar(categorical({'T1_1', 'T1_2', 'T1_3', 'T1_4', 'T2_1', 'T2_2', 'T2_3', 'T2_4'},...
        {'T1_1', 'T2_1', 'T1_2', 'T2_2', 'T1_3', 'T2_3', 'T1_4', 'T2_4'}),meanAll);


    mysave(f, fullfile(outputfolder{1}, 'R2MaxSmallVsBig.txt'))
    [p1, h1] = ranksum(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation1), tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation2));
    [p2, h2] = ranksum(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation1), tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation2));
    [p3, h3] = ranksum(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation1), tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation2));
    [p4, h4] = ranksum(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation1), tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation2));

    text = '';
    text = strcat(text, sprintf('Cluster1 h1 %f, p %f \\n', h1,p1));
    text = strcat(text, sprintf('Cluster2 h1 %f, p %f \\n', h2,p2));
    text = strcat(text, sprintf('Cluster3 h1 %f, p %f \\n', h3,p3));
    text = strcat(text, sprintf('Cluster4 h1 %f, p %f \\n', h4,p4));

    fid=fopen(fullfile(outputfolder{1}, 'statisticranksumR2MaxHemiSmallVsBig.txt'),'w');
    fprintf(fid, text);
    fclose(fid);
    
    calcAnovaForPearsonCompare(matAll, outputfolder{1});
    calcAnovaForPearsonCompareALL(matAll, outputfolder{1})
%     calcAnovaForPearsonCompare2(matAll2, outputfolder{1});
end

function calcAnovaForPearsonCompare(matAll, outputfolder)
    for i = 1:4
        valuesIndex = find(matAll(:, 2) == i);
        [p,tbl,stats] = anovan(matAll(valuesIndex, 1),{matAll(valuesIndex, 4) matAll(valuesIndex, 3)},'model','interaction','varnames',{'B/W','Type'});
        writetable(cell2table(tbl), fullfile(outputfolder, sprintf('AnovaTestByCluster_%d_TOPearsonValuesB_W.csv', i)));
        f = figure; hold on;
        
        y1 = matAll(matAll(:, 2) == i & matAll(:, 3) == 1& matAll(:, 4) == 1, 1);
        y2 = matAll(matAll(:, 2) == i & matAll(:, 3) == 1& matAll(:, 4) == 2, 1);
        y3 = matAll(matAll(:, 2) == i & matAll(:, 3) == 2& matAll(:, 4) == 1, 1);
        y4 = matAll(matAll(:, 2) == i & matAll(:, 3) == 2& matAll(:, 4) == 2, 1);
        meanAll1 = [mean(y1), mean(y2)];
        meanAll2 = [mean(y3), mean(y4)];
        
        c1 = cell(1, length(y1));
        c1(:) = {'Big-Within'};
        c2 = cell(1, length(y2));
        c2(:) = {'Big-Between'};
        c3 = cell(1, length(y3));
        c3(:) = {'Small-Within'};
        c4 = cell(1, length(y4));
        c4(:) = {'Small-Between'};
        
        bc1 = boxchart(categorical(c1),y1);
        bc2 = boxchart(categorical(c2),y2);
        bc3 = boxchart(categorical(c3),y3);
        bc4 = boxchart(categorical(c4),y4);
        
        plot(categorical({'Big-Within', 'Big-Between'}), meanAll1, '-*k');
        plot(categorical({'Small-Within', 'Small-Between'}), meanAll2, '-*k');
        
        mysave(f, fullfile(outputfolder, sprintf('AnovaTestByCluster_%d_TOPearsonValuesB_WFig', i)));
        
        textAnova = '';
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Big = %f, std %f \\n ',i, mean(y1),  std(y1)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Big = %f, std %f \\n ',i, mean(y2), std(y2)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Small = %f, std %f \\n ',i, mean(y3), std(y3)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Small = %f, std %f \\n ',i, mean(y4), std(y4)));
        
        fid=fopen(fullfile(outputfolder, sprintf('Mean_PearsonWB_cluster%d.txt', i)),'w');
        fprintf(fid, textAnova);
        fclose(fid);
        
        textAnova = '';
        [c,~,~,groupnames] = multcompare(stats,'Display', 'off', 'Dimension', [1,2]);
        for j = 1:size(c, 1)
            textAnova = strcat(textAnova, sprintf('%s vs %s - pValue : %.4f, mean Diff : %.4f, CI_L : %.4f, CI_H : %.4f \\n ', ...
                groupnames{c(j, 1)}, groupnames{c(j, 2)}, c(j, 6), c(j, 4), c(j, 3), c(j, 5)));
        end

        fid=fopen(fullfile(outputfolder, sprintf('CompareGroups_PearsonWB_cluster%d.txt', i)),'w');
        fprintf(fid, textAnova);
        fclose(fid);
    end
end

function calcAnovaForPearsonCompareALL(matAll, outputfolder)
    
    valuesIndex = 1:length(matAll(:, 2));
    [p,tbl,stats] = anovan(matAll(valuesIndex, 1),{matAll(valuesIndex, 4) matAll(valuesIndex, 3)},'model','interaction','varnames',{'B/W','Type'});
    writetable(cell2table(tbl), fullfile(outputfolder, sprintf('AnovaTestByCluster_%d_TOPearsonValuesB_W.csv', 0)));
    f = figure; hold on;

    y1 = matAll(matAll(:, 3) == 1& matAll(:, 4) == 1, 1);
    y2 = matAll(matAll(:, 3) == 1& matAll(:, 4) == 2, 1);
    y3 = matAll(matAll(:, 3) == 2& matAll(:, 4) == 1, 1);
    y4 = matAll(matAll(:, 3) == 2& matAll(:, 4) == 2, 1);
    meanAll1 = [mean(y1), mean(y2)];
    meanAll2 = [mean(y3), mean(y4)];

    c1 = cell(1, length(y1));
    c1(:) = {'Big-Within'};
    c2 = cell(1, length(y2));
    c2(:) = {'Big-Between'};
    c3 = cell(1, length(y3));
    c3(:) = {'Small-Within'};
    c4 = cell(1, length(y4));
    c4(:) = {'Small-Between'};

    bc1 = boxchart(categorical(c1),y1);
    bc2 = boxchart(categorical(c2),y2);
    bc3 = boxchart(categorical(c3),y3);
    bc4 = boxchart(categorical(c4),y4);

    plot(categorical({'Big-Within', 'Big-Between'}), meanAll1, '-*k');
    plot(categorical({'Small-Within', 'Small-Between'}), meanAll2, '-*k');

    mysave(f, fullfile(outputfolder, sprintf('AnovaTestByCluster_%d_TOPearsonValuesB_WFig', 0)));

    textAnova = '';
    textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Big = %f, std %f \\n ',i, mean(y1),  std(y1)));
    textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Big = %f, std %f \\n ',i, mean(y2), std(y2)));
    textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Small = %f, std %f \\n ',i, mean(y3), std(y3)));
    textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Small = %f, std %f \\n ',i, mean(y4), std(y4)));

    fid=fopen(fullfile(outputfolder, sprintf('Mean_PearsonWB_cluster%d.txt', 0)),'w');
    fprintf(fid, textAnova);
    fclose(fid);

    textAnova = '';
    [c,~,~,groupnames] = multcompare(stats,'Display', 'off', 'Dimension', [1,2]);
    for j = 1:size(c, 1)
        textAnova = strcat(textAnova, sprintf('%s vs %s - pValue : %.4f, mean Diff : %.4f, CI_L : %.4f, CI_H : %.4f \\n ', ...
            groupnames{c(j, 1)}, groupnames{c(j, 2)}, c(j, 6), c(j, 4), c(j, 3), c(j, 5)));
    end

    fid=fopen(fullfile(outputfolder, sprintf('CompareGroups_PearsonWB_cluster%d.txt', 0)),'w');
    fprintf(fid, textAnova);
    fclose(fid);

end


function calcAnovaForPearsonCompare2(matAll, outputfolder)
    for i = 1:4
        valuesIndex = find(matAll(:, 2) == i);
        [p,tbl,stats] = anovan(matAll(valuesIndex, 1),{matAll(valuesIndex, 4) matAll(valuesIndex, 3)},'model','interaction','varnames',{'BSec/BMain/W','Type'});
        writetable(cell2table(tbl), fullfile(outputfolder, sprintf('AnovaTestByCluster_%d_TOPearsonValuesB_W_mainvsSec.csv', i)));
        f = figure; hold on;
        
        y1 = matAll(matAll(:, 2) == i & matAll(:, 3) == 1& matAll(:, 4) == 1, 1);
        y2 = matAll(matAll(:, 2) == i & matAll(:, 3) == 1& matAll(:, 4) == 2, 1);
        y5 = matAll(matAll(:, 2) == i & matAll(:, 3) == 1& matAll(:, 4) == 3, 1);
        y3 = matAll(matAll(:, 2) == i & matAll(:, 3) == 2& matAll(:, 4) == 1, 1);
        y4 = matAll(matAll(:, 2) == i & matAll(:, 3) == 2& matAll(:, 4) == 2, 1);
        y6 = matAll(matAll(:, 2) == i & matAll(:, 3) == 2& matAll(:, 4) == 3, 1);
        meanAll1 = [mean(y1), mean(y5),mean(y2)];
        meanAll2 = [mean(y3), mean(y6),mean(y4)];
        
        c1 = cell(1, length(y1));
        c1(:) = {'Big-Within'};
        c2 = cell(1, length(y2));
        c2(:) = {'Big-Between-Main'};
        c5 = cell(1, length(y5));
        c5(:) = {'Big-Between-Sec'};
        c3 = cell(1, length(y3));
        c3(:) = {'Small-Within'};
        c4 = cell(1, length(y4));
        c4(:) = {'Small-Between-Main'};
        c6 = cell(1, length(y6));
        c6(:) = {'Small-Between-Sec'};
        
        bc1 = boxchart(categorical(c1),y1);
        bc5 = boxchart(categorical(c5),y5);
        bc2 = boxchart(categorical(c2),y2);
        bc3 = boxchart(categorical(c3),y3);
        bc6 = boxchart(categorical(c6),y6);
        bc4 = boxchart(categorical(c4),y4);
        
        plot(categorical({'Big-Within', 'Big-Between-Sec', 'Big-Between-Main'}), meanAll1, '-*k');
        plot(categorical({'Small-Within', 'Small-Between-Sec', 'Small-Between-Main'}), meanAll2, '-*k');
        
        mysave(f, fullfile(outputfolder, sprintf('AnovaTestByCluster_%d_TOPearsonValuesB_W_secvsmainFig', i)));
        
        textAnova = '';
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Big = %f \\n ',i, mean(y1)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Main Big = %f \\n ',i, mean(y2)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Sec Big = %f \\n ',i, mean(y5)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Small = %f \\n ',i, mean(y3)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Main Small = %f \\n ',i, mean(y4)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between sec Small = %f \\n ',i, mean(y6)));
        
        fid=fopen(fullfile(outputfolder, sprintf('Mean_PearsonWB_secMain_cluster%d.txt', i)),'w');
        fprintf(fid, textAnova);
        fclose(fid);
        
        textAnova = '';
        [c,~,~,groupnames] = multcompare(stats,'Display', 'off', 'Dimension', [1,2]);
        for j = 1:size(c, 1)
            textAnova = strcat(textAnova, sprintf('%s vs %s - pValue : %.4f, mean Diff : %.4f, CI_L : %.4f, CI_H : %.4f \\n ', ...
                groupnames{c(j, 1)}, groupnames{c(j, 2)}, c(j, 6), c(j, 4), c(j, 3), c(j, 5)));
        end

        fid=fopen(fullfile(outputfolder, sprintf('CompareGroups_PearsonWB_secMain_cluster%d.txt', i)),'w');
        fprintf(fid, textAnova);
        fclose(fid);
    end
end

function [matAll, bslop] = BetweenVsWithinR2(indexLocation,indexHemi, tableResults, outputfolder, alphValue)
    cluster1betweenAll = [];
    cluster1withinAll = [];
    cluster2betweenAll = [];
    cluster2withinAll = [];
    cluster3betweenAll = [];
    cluster3withinAll = [];
    cluster4betweenAll = [];
    cluster4withinAll = [];

    tableResults.b1_cluster1 = nan(size(tableResults, 1),1);
    tableResults.b1_cluster2 = nan(size(tableResults, 1),1);
    tableResults.b1_cluster3 = nan(size(tableResults, 1),1);
    tableResults.b1_cluster4 = nan(size(tableResults, 1),1);
    
    for i = 1:size(indexLocation, 1)
        cluster1File = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*.csv'];
        
        fileList = dir(cluster1File);
        
        xVals = [];
        yVals = [];
        x1Vals = {};
        y1Vals = {};
        in_k = 1;
        
        temp1 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
            
            tableR.Var2 = power(tableR.Var2, 2);
            
            if contains(fileList(j).name, 'BetweenMainDepth')
                cluster1betweenAll(end+1) = nanmean(tableR.Var2);
            else
                temp1(end+1:end+length(tableR.Var2)) = tableR.Var2;
                x1Vals(in_k) = {tableR.Var1};
                y1Vals(in_k) = {tableR.Var2};
                 in_k = in_k + 1;
            end
            
            
            xVals(end+1:end+length(tableR.Var1)) = tableR.Var1;
            yVals(end+1:end+length(tableR.Var2)) = tableR.Var2;
        end
        cluster1withinAll(end+1) = nanmean(temp1);
            
        cluster1withinAll(isnan(cluster1withinAll)) = [];
        cluster1betweenAll(isnan(cluster1betweenAll)) = [];
         
        mdAll1 = fitglm(xVals ./ max(xVals), yVals);
        b1_cluster1 = mdAll1.Coefficients.Estimate(2);
        
        tableResults.b1_cluster1(indexLocation(i)) = b1_cluster1;
        
        mdSub1 = fitglm(x1Vals{1} ./ max(x1Vals{1}), y1Vals{1});
        b1_sec1_cluster1 = mdSub1.Coefficients.Estimate(2);
        tableResults.b1_sub1_cluster1(indexLocation(i)) = b1_sec1_cluster1;
        
        mdSub2 = fitglm(x1Vals{2} ./ max(x1Vals{2}), y1Vals{2});
        b1_sec2_cluster1 = mdSub2.Coefficients.Estimate(2);
        tableResults.b1_sub2_cluster1(indexLocation(i)) = b1_sec2_cluster1;
        
        cluster2File = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*.csv'];
        
        fileList = dir(cluster2File);
        
        xVals = [];
        yVals = [];
       
        x1Vals = {};
        y1Vals = {};
        in_k = 1;
        
        temp2 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
            
            tableR.Var2 = power(tableR.Var2, 2);
            if contains(fileList(j).name, 'BetweenMainDepth')
                cluster2betweenAll(end+1) = nanmean(tableR.Var2);
            else
                temp2(end+1:end+length(tableR.Var2)) = tableR.Var2;
                x1Vals(in_k) = {tableR.Var1};
                y1Vals(in_k) = {tableR.Var2};
                in_k = in_k + 1;
            
            end
            
            
            xVals(end+1:end+length(tableR.Var1)) = tableR.Var1;
            yVals(end+1:end+length(tableR.Var2)) = tableR.Var2;
        end
        cluster2withinAll(end+1) = nanmean(temp2);
            
        cluster2withinAll(isnan(cluster2withinAll)) = [];
        cluster2betweenAll(isnan(cluster2betweenAll)) = [];
         
        
        mdAll1 = fitglm(xVals ./ max(xVals), yVals);
        b1_cluster2 = mdAll1.Coefficients.Estimate(2);
        
        tableResults.b1_cluster2(indexLocation(i)) = b1_cluster2;
        
        
        mdSub1 = fitglm(x1Vals{1} ./ max(x1Vals{1}), y1Vals{1});
        b1_sec1_cluster2 = mdSub1.Coefficients.Estimate(2);
        tableResults.b1_sub1_cluster2(indexLocation(i)) = b1_sec1_cluster2;
        
        mdSub2 = fitglm(x1Vals{2} ./ max(x1Vals{2}), y1Vals{2});
        b1_sec2_cluster2 = mdSub2.Coefficients.Estimate(2);
        tableResults.b1_sub2_cluster2(indexLocation(i)) = b1_sec2_cluster2;
        
        
        cluster3File = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*.csv'];
        
        fileList = dir(cluster3File);
        
        xVals = [];
        yVals = [];
       
        x1Vals = {};
        y1Vals = {};
        in_k = 1;
      
        temp3 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
            
            tableR.Var2 = power(tableR.Var2, 2);
            if contains(fileList(j).name, 'BetweenMainDepth')
                cluster3betweenAll(end+1) = nanmean(tableR.Var2);
            else
                temp3(end+1:end+length(tableR.Var2)) = tableR.Var2;
                x1Vals(in_k) = {tableR.Var1};
                y1Vals(in_k) = {tableR.Var2};
                in_k = in_k + 1;
            
            end
            
            
            xVals(end+1:end+length(tableR.Var1)) = tableR.Var1;
            yVals(end+1:end+length(tableR.Var2)) = tableR.Var2;
        end
        cluster3withinAll(end+1) = nanmean(temp3);
            
        cluster3withinAll(isnan(cluster3withinAll)) = [];
        cluster3betweenAll(isnan(cluster3betweenAll)) = [];
           
        
        mdAll1 = fitglm(xVals ./ max(xVals), yVals);
        b1_cluster3 = mdAll1.Coefficients.Estimate(2);
        
        tableResults.b1_cluster3(indexLocation(i)) = b1_cluster3;
        
        
        mdSub1 = fitglm(x1Vals{1} ./ max(x1Vals{1}), y1Vals{1});
        b1_sec1_cluster3 = mdSub1.Coefficients.Estimate(2);
        tableResults.b1_sub1_cluster3(indexLocation(i)) = b1_sec1_cluster3;
        
        mdSub2 = fitglm(x1Vals{2} ./ max(x1Vals{2}), y1Vals{2});
        b1_sec2_cluster3 = mdSub2.Coefficients.Estimate(2);
        tableResults.b1_sub2_cluster3(indexLocation(i)) = b1_sec2_cluster3;
        
        
        cluster4File = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*.csv'];
        
        fileList = dir(cluster4File);
        
        xVals = [];
        yVals = [];
       
        x1Vals = {};
        y1Vals = {};
        in_k = 1;
      
        temp4 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
            
            tableR.Var2 = power(tableR.Var2, 2);
            if contains(fileList(j).name, 'BetweenMainDepth')
                cluster4betweenAll(end+1) = nanmean(tableR.Var2);
            else
                temp4(end+1:end+length(tableR.Var2)) = tableR.Var2;
                x1Vals(in_k) = {tableR.Var1};
                y1Vals(in_k) = {tableR.Var2};
                in_k = in_k + 1;
            
            end
            
            
            xVals(end+1:end+length(tableR.Var1)) = tableR.Var1;
            yVals(end+1:end+length(tableR.Var2)) = tableR.Var2;
        end
        
        cluster4withinAll(end+1) = nanmean(temp4);
        
        cluster4withinAll(isnan(cluster4withinAll)) = [];
        cluster4betweenAll(isnan(cluster4betweenAll)) = [];
        
        mdAll1 = fitglm(xVals ./ max(xVals), yVals);
        b1_cluster4 = mdAll1.Coefficients.Estimate(2);
        
        tableResults.b1_cluster4(indexLocation(i)) = b1_cluster4;
        
        mdSub1 = fitglm(x1Vals{1} ./ max(x1Vals{1}), y1Vals{1});
        b1_sec1_cluster1 = mdSub1.Coefficients.Estimate(2);
        tableResults.b1_sub1_cluster4(indexLocation(i)) = b1_sec1_cluster1;
        
        mdSub2 = fitglm(x1Vals{2} ./ max(x1Vals{2}), y1Vals{2});
        b1_sec2_cluster4 = mdSub2.Coefficients.Estimate(2);
        tableResults.b1_sub2_cluster4(indexLocation(i)) = b1_sec2_cluster4;
        
        
    end
    
    save(fullfile(outputfolder, 'tablewithB1'), 'tableResults')
    
    alphaFix1 = alphValue / (sqrt(length(cluster1withinAll) / 100));
    alphaFix2 = alphValue / (sqrt(length(cluster2withinAll) / 100));
    alphaFix3 = alphValue / (sqrt(length(cluster3withinAll) / 100));
    alphaFix4 = alphValue / (sqrt(length(cluster4withinAll) / 100));
    
    [h1,p1] = ttest2(cluster1withinAll, cluster1betweenAll, 'Alpha', alphaFix1);
    [h2,p2] = ttest2(cluster2withinAll, cluster2betweenAll, 'Alpha', alphaFix2);    
    [h3,p3] = ttest2(cluster3withinAll, cluster3betweenAll, 'Alpha', alphaFix3);
    [h4,p4] = ttest2(cluster4withinAll, cluster4betweenAll, 'Alpha', alphaFix4);
       
    tableTtest = table([h1;h2;h3;h4], [p1;p2;p3;p4], 'VariableNames', {'H', 'Pvalue'}, 'RowNames', {'1','2','3','4'});
    writetable(tableTtest, fullfile(outputfolder, 'TTestResultsWithinBetweenByData.csv'));
   
    f = figure; hold on;
    groupsBox1 = {};
    groupsBox1(end+1:end+length(cluster1withinAll)) = {'Within-Cluster1'};
    groupsBox1(end+1:end+length(cluster2withinAll)) = {'Within-Cluster2'};
    groupsBox1(end+1:end+length(cluster3withinAll)) = {'Within-Cluster3'};
    groupsBox1(end+1:end+length(cluster4withinAll)) = {'Within-Cluster4'};
    groupsBox2 = {};
    groupsBox2(end+1:end+length(cluster1betweenAll)) = {'Between-Cluster1'};
    groupsBox2(end+1:end+length(cluster2betweenAll)) = {'Between-Cluster2'};
    groupsBox2(end+1:end+length(cluster3betweenAll)) = {'Between-Cluster3'};
    groupsBox2(end+1:end+length(cluster4betweenAll)) = {'Between-Cluster4'};
    
    meanAll(1) = mean(cluster1withinAll);
    meanAll(2) = mean(cluster2withinAll);
    meanAll(3) = mean(cluster3withinAll);   
    meanAll(4) = mean(cluster4withinAll);  
    meanAll(5) = mean(cluster1betweenAll);
    meanAll(6) = mean(cluster2betweenAll);    
    meanAll(7) = mean(cluster3betweenAll);  
    meanAll(8) = mean(cluster4betweenAll);
    

    bC = boxchart(categorical(groupsBox1),[cluster1withinAll'; ...
        cluster2withinAll';...
        cluster3withinAll';...
        cluster4withinAll']);
    
    bC2 = boxchart(categorical(groupsBox2),[cluster1betweenAll';...
        cluster2betweenAll';...
        cluster3betweenAll';...
        cluster4betweenAll']);
    
  
    bC.BoxFaceColor = [0,0,0];
    bC.BoxFaceAlpha = 0.4;
    bC.MarkerColor = [0,0,0];
    bC2.BoxFaceColor = [0,0,255] ./ 255;
    bC2.MarkerColor = [0,0,255] ./ 255;
    bC2.BoxFaceAlpha = 0.4;
    
    plot(meanAll, '-*k');

    mysave(f, fullfile(outputfolder, 'PlotMeanAndSTDByDataWithinBetweenAllClusters'));
    
    matResultsV = [cluster1withinAll';cluster1betweenAll';cluster2withinAll';cluster2betweenAll';...
        cluster3withinAll';cluster3betweenAll';cluster4withinAll';cluster4betweenAll'];
    
    matResultsCluster = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenAll),1);...
        ones(length(cluster2withinAll),1)*2;ones(length(cluster2betweenAll),1)*2;...
        ones(length(cluster3withinAll),1)*3;ones(length(cluster3betweenAll),1)*3;...
        ones(length(cluster4withinAll),1)*4; ones(length(cluster4betweenAll),1)*4];
    
    matResultsType = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenAll),1);...
        ones(length(cluster2withinAll),1);ones(length(cluster2betweenAll),1);...
        ones(length(cluster3withinAll),1);ones(length(cluster3betweenAll),1);...
        ones(length(cluster4withinAll),1); ones(length(cluster4betweenAll),1)];
    
    matResultsBW = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenAll),1)*2;...
        ones(length(cluster2withinAll),1);ones(length(cluster2betweenAll),1)*2;...
        ones(length(cluster3withinAll),1);ones(length(cluster3betweenAll),1)*2;...
        ones(length(cluster4withinAll),1); ones(length(cluster4betweenAll),1)*2];
    
    matAll = [matResultsV, matResultsCluster, matResultsType, matResultsBW];
    
    f = figure; hold on;
    title('Linear regression slop (b1)');
    groupsBox1 = {};
    groupsBox1(end+1:end+length(indexLocation)) = {'Cluster1'};
    groupsBox1(end+1:end+length(indexLocation)) = {'Cluster2'};
    groupsBox1(end+1:end+length(indexLocation)) = {'Cluster3'};
    groupsBox1(end+1:end+length(indexLocation)) = {'Cluster4'};
    meanAll = [];
    bC = boxchart(categorical(groupsBox1), [tableResults.b1_cluster1(indexLocation); ...
        tableResults.b1_cluster2(indexLocation);...
        tableResults.b1_cluster3(indexLocation);...
        tableResults.b1_cluster4(indexLocation)]);
    meanAll(1) = mean(tableResults.b1_cluster1(indexLocation));
    meanAll(2) = mean(tableResults.b1_cluster2(indexLocation));
    meanAll(3) = mean(tableResults.b1_cluster3(indexLocation));   
    meanAll(4) = mean(tableResults.b1_cluster4(indexLocation));  
    plot(meanAll, '-*k');
    ylim([-1.2, 0.2]);
    mysave(f, fullfile(outputfolder, 'B1_slop'));   
    
    f = figure; hold on;
    title('Linear regression slop (b1) Hemi');
    groupsBox1 = {};
    groupsBox1(end+1:end+length(indexHemi)) = {'Cluster1'};
    groupsBox1(end+1:end+length(indexHemi)) = {'Cluster2'};
    groupsBox1(end+1:end+length(indexHemi)) = {'Cluster3'};
    groupsBox1(end+1:end+length(indexHemi)) = {'Cluster4'};
    meanAll = [];
    bC = boxchart(categorical(groupsBox1),...
        [min([tableResults.b1_sub1_cluster1(indexHemi),tableResults.b1_sub2_cluster1(indexHemi)],[],2); ...
        min([tableResults.b1_sub1_cluster2(indexHemi),tableResults.b1_sub2_cluster2(indexHemi)],[], 2);...
        min([tableResults.b1_sub1_cluster3(indexHemi),tableResults.b1_sub2_cluster3(indexHemi)],[], 2);...
        min([tableResults.b1_sub1_cluster4(indexHemi),tableResults.b1_sub2_cluster4(indexHemi)],[], 2)]);
    meanAll(1) = mean(min([tableResults.b1_sub1_cluster1(indexHemi),tableResults.b1_sub2_cluster1(indexHemi)],[],2));
    meanAll(2) = mean(min([tableResults.b1_sub1_cluster2(indexHemi),tableResults.b1_sub2_cluster2(indexHemi)],[],2));
    meanAll(3) = mean(min([tableResults.b1_sub1_cluster3(indexHemi),tableResults.b1_sub2_cluster3(indexHemi)],[],2));   
    meanAll(4) = mean(min([tableResults.b1_sub1_cluster4(indexHemi),tableResults.b1_sub2_cluster4(indexHemi)],[], 2));  
    plot(meanAll, '-*k');
    ylim([-1.2, 0.2]);
    mysave(f, fullfile(outputfolder, 'B1_slop_Hemi'));  
    
    statPlot([tableResults.b1_cluster1(indexLocation), ...
        tableResults.b1_cluster2(indexLocation),...
        tableResults.b1_cluster3(indexLocation),...
        tableResults.b1_cluster4(indexLocation)], 'slop(b1)', 'Linear regression slop (b1)', outputfolder, [0,0,128]./255);
    
    bslop = [tableResults.b1_cluster1(indexLocation), ...
        tableResults.b1_cluster2(indexLocation),...
        tableResults.b1_cluster3(indexLocation),...
        tableResults.b1_cluster4(indexLocation)];
end

function [matAll, bslop] = BetweenVsWithin(indexLocation,indexHemi, tableResults, outputfolder, alphValue)
    cluster1betweenAll = [];
    cluster1withinAll = [];
    cluster2betweenAll = [];
    cluster2withinAll = [];
    cluster3betweenAll = [];
    cluster3withinAll = [];
    cluster4betweenAll = [];
    cluster4withinAll = [];

    tableResults.b1_cluster1 = nan(size(tableResults, 1),1);
    tableResults.b1_cluster2 = nan(size(tableResults, 1),1);
    tableResults.b1_cluster3 = nan(size(tableResults, 1),1);
    tableResults.b1_cluster4 = nan(size(tableResults, 1),1);
    
    for i = 1:size(indexLocation, 1)
        cluster1File = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*.csv'];
        
        fileList = dir(cluster1File);
        
        xVals = [];
        yVals = [];
        x1Vals = {};
        y1Vals = {};
        in_k = 1;
        
        temp1 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
            
            if contains(fileList(j).name, 'BetweenMainDepth')
                cluster1betweenAll(end+1) = nanmean(tableR.Var2);
            else
                temp1(end+1:end+length(tableR.Var2)) = tableR.Var2;
                x1Vals(in_k) = {tableR.Var1};
                y1Vals(in_k) = {tableR.Var2};
                 in_k = in_k + 1;
            end
            
            
            xVals(end+1:end+length(tableR.Var1)) = tableR.Var1;
            yVals(end+1:end+length(tableR.Var2)) = tableR.Var2;
        end
        cluster1withinAll(end+1) = nanmean(temp1);
            
        cluster1withinAll(isnan(cluster1withinAll)) = [];
        cluster1betweenAll(isnan(cluster1betweenAll)) = [];
         
        mdAll1 = fitglm(xVals ./ max(xVals), yVals);
        b1_cluster1 = mdAll1.Coefficients.Estimate(2);
        
        tableResults.b1_cluster1(indexLocation(i)) = b1_cluster1;
        
        mdSub1 = fitglm(x1Vals{1} ./ max(x1Vals{1}), y1Vals{1});
        b1_sec1_cluster1 = mdSub1.Coefficients.Estimate(2);
        tableResults.b1_sub1_cluster1(indexLocation(i)) = b1_sec1_cluster1;
        
        mdSub2 = fitglm(x1Vals{2} ./ max(x1Vals{2}), y1Vals{2});
        b1_sec2_cluster1 = mdSub2.Coefficients.Estimate(2);
        tableResults.b1_sub2_cluster1(indexLocation(i)) = b1_sec2_cluster1;
        
        cluster2File = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*.csv'];
        
        fileList = dir(cluster2File);
        
        xVals = [];
        yVals = [];
       
        x1Vals = {};
        y1Vals = {};
        in_k = 1;
        
        temp2 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
            
            if contains(fileList(j).name, 'BetweenMainDepth')
                cluster2betweenAll(end+1) = nanmean(tableR.Var2);
            else
                temp2(end+1:end+length(tableR.Var2)) = tableR.Var2;
                x1Vals(in_k) = {tableR.Var1};
                y1Vals(in_k) = {tableR.Var2};
                in_k = in_k + 1;
            
            end
            
            
            xVals(end+1:end+length(tableR.Var1)) = tableR.Var1;
            yVals(end+1:end+length(tableR.Var2)) = tableR.Var2;
        end
        cluster2withinAll(end+1) = nanmean(temp2);
            
        cluster2withinAll(isnan(cluster2withinAll)) = [];
        cluster2betweenAll(isnan(cluster2betweenAll)) = [];
         
        
        mdAll1 = fitglm(xVals ./ max(xVals), yVals);
        b1_cluster2 = mdAll1.Coefficients.Estimate(2);
        
        tableResults.b1_cluster2(indexLocation(i)) = b1_cluster2;
        
        
        mdSub1 = fitglm(x1Vals{1} ./ max(x1Vals{1}), y1Vals{1});
        b1_sec1_cluster2 = mdSub1.Coefficients.Estimate(2);
        tableResults.b1_sub1_cluster2(indexLocation(i)) = b1_sec1_cluster2;
        
        mdSub2 = fitglm(x1Vals{2} ./ max(x1Vals{2}), y1Vals{2});
        b1_sec2_cluster2 = mdSub2.Coefficients.Estimate(2);
        tableResults.b1_sub2_cluster2(indexLocation(i)) = b1_sec2_cluster2;
        
        
        cluster3File = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*.csv'];
        
        fileList = dir(cluster3File);
        
        xVals = [];
        yVals = [];
       
        x1Vals = {};
        y1Vals = {};
        in_k = 1;
      
        temp3 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
            
            if contains(fileList(j).name, 'BetweenMainDepth')
                cluster3betweenAll(end+1) = nanmean(tableR.Var2);
            else
                temp3(end+1:end+length(tableR.Var2)) = tableR.Var2;
                x1Vals(in_k) = {tableR.Var1};
                y1Vals(in_k) = {tableR.Var2};
                in_k = in_k + 1;
            
            end
            
            
            xVals(end+1:end+length(tableR.Var1)) = tableR.Var1;
            yVals(end+1:end+length(tableR.Var2)) = tableR.Var2;
        end
        cluster3withinAll(end+1) = nanmean(temp3);
            
        cluster3withinAll(isnan(cluster3withinAll)) = [];
        cluster3betweenAll(isnan(cluster3betweenAll)) = [];
           
        
        mdAll1 = fitglm(xVals ./ max(xVals), yVals);
        b1_cluster3 = mdAll1.Coefficients.Estimate(2);
        
        tableResults.b1_cluster3(indexLocation(i)) = b1_cluster3;
        
        
        mdSub1 = fitglm(x1Vals{1} ./ max(x1Vals{1}), y1Vals{1});
        b1_sec1_cluster3 = mdSub1.Coefficients.Estimate(2);
        tableResults.b1_sub1_cluster3(indexLocation(i)) = b1_sec1_cluster3;
        
        mdSub2 = fitglm(x1Vals{2} ./ max(x1Vals{2}), y1Vals{2});
        b1_sec2_cluster3 = mdSub2.Coefficients.Estimate(2);
        tableResults.b1_sub2_cluster3(indexLocation(i)) = b1_sec2_cluster3;
        
        
        cluster4File = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*.csv'];
        
        fileList = dir(cluster4File);
        
        xVals = [];
        yVals = [];
       
        x1Vals = {};
        y1Vals = {};
        in_k = 1;
      
        temp4 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
            
            if contains(fileList(j).name, 'BetweenMainDepth')
                cluster4betweenAll(end+1) = nanmean(tableR.Var2);
            else
                temp4(end+1:end+length(tableR.Var2)) = tableR.Var2;
                x1Vals(in_k) = {tableR.Var1};
                y1Vals(in_k) = {tableR.Var2};
                in_k = in_k + 1;
            
            end
            
            
            xVals(end+1:end+length(tableR.Var1)) = tableR.Var1;
            yVals(end+1:end+length(tableR.Var2)) = tableR.Var2;
        end
        
        cluster4withinAll(end+1) = nanmean(temp4);
        
        cluster4withinAll(isnan(cluster4withinAll)) = [];
        cluster4betweenAll(isnan(cluster4betweenAll)) = [];
        
        mdAll1 = fitglm(xVals ./ max(xVals), yVals);
        b1_cluster4 = mdAll1.Coefficients.Estimate(2);
        
        tableResults.b1_cluster4(indexLocation(i)) = b1_cluster4;
        
        mdSub1 = fitglm(x1Vals{1} ./ max(x1Vals{1}), y1Vals{1});
        b1_sec1_cluster1 = mdSub1.Coefficients.Estimate(2);
        tableResults.b1_sub1_cluster4(indexLocation(i)) = b1_sec1_cluster1;
        
        mdSub2 = fitglm(x1Vals{2} ./ max(x1Vals{2}), y1Vals{2});
        b1_sec2_cluster4 = mdSub2.Coefficients.Estimate(2);
        tableResults.b1_sub2_cluster4(indexLocation(i)) = b1_sec2_cluster4;
        
        
    end
    
    save(fullfile(outputfolder, 'tablewithB1'), 'tableResults')
    
    alphaFix1 = alphValue / (sqrt(length(cluster1withinAll) / 100));
    alphaFix2 = alphValue / (sqrt(length(cluster2withinAll) / 100));
    alphaFix3 = alphValue / (sqrt(length(cluster3withinAll) / 100));
    alphaFix4 = alphValue / (sqrt(length(cluster4withinAll) / 100));
    
    [h1,p1] = ttest2(cluster1withinAll, cluster1betweenAll, 'Alpha', alphaFix1);
    [h2,p2] = ttest2(cluster2withinAll, cluster2betweenAll, 'Alpha', alphaFix2);    
    [h3,p3] = ttest2(cluster3withinAll, cluster3betweenAll, 'Alpha', alphaFix3);
    [h4,p4] = ttest2(cluster4withinAll, cluster4betweenAll, 'Alpha', alphaFix4);
       
    tableTtest = table([h1;h2;h3;h4], [p1;p2;p3;p4], 'VariableNames', {'H', 'Pvalue'}, 'RowNames', {'1','2','3','4'});
    writetable(tableTtest, fullfile(outputfolder, 'TTestResultsWithinBetweenByData.csv'));
   
    f = figure; hold on;
    groupsBox1 = {};
    groupsBox1(end+1:end+length(cluster1withinAll)) = {'Within-Cluster1'};
    groupsBox1(end+1:end+length(cluster2withinAll)) = {'Within-Cluster2'};
    groupsBox1(end+1:end+length(cluster3withinAll)) = {'Within-Cluster3'};
    groupsBox1(end+1:end+length(cluster4withinAll)) = {'Within-Cluster4'};
    groupsBox2 = {};
    groupsBox2(end+1:end+length(cluster1betweenAll)) = {'Between-Cluster1'};
    groupsBox2(end+1:end+length(cluster2betweenAll)) = {'Between-Cluster2'};
    groupsBox2(end+1:end+length(cluster3betweenAll)) = {'Between-Cluster3'};
    groupsBox2(end+1:end+length(cluster4betweenAll)) = {'Between-Cluster4'};
    
    meanAll(1) = mean(cluster1withinAll);
    meanAll(2) = mean(cluster2withinAll);
    meanAll(3) = mean(cluster3withinAll);   
    meanAll(4) = mean(cluster4withinAll);  
    meanAll(5) = mean(cluster1betweenAll);
    meanAll(6) = mean(cluster2betweenAll);    
    meanAll(7) = mean(cluster3betweenAll);  
    meanAll(8) = mean(cluster4betweenAll);
    

    bC = boxchart(categorical(groupsBox1),[cluster1withinAll'; ...
        cluster2withinAll';...
        cluster3withinAll';...
        cluster4withinAll']);
    
    bC2 = boxchart(categorical(groupsBox2),[cluster1betweenAll';...
        cluster2betweenAll';...
        cluster3betweenAll';...
        cluster4betweenAll']);
    
  
    bC.BoxFaceColor = [0,0,0];
    bC.BoxFaceAlpha = 0.4;
    bC.MarkerColor = [0,0,0];
    bC2.BoxFaceColor = [0,0,255] ./ 255;
    bC2.MarkerColor = [0,0,255] ./ 255;
    bC2.BoxFaceAlpha = 0.4;
    
    plot(meanAll, '-*k');

    mysave(f, fullfile(outputfolder, 'PlotMeanAndSTDByDataWithinBetweenAllClusters'));
    
    matResultsV = [cluster1withinAll';cluster1betweenAll';cluster2withinAll';cluster2betweenAll';...
        cluster3withinAll';cluster3betweenAll';cluster4withinAll';cluster4betweenAll'];
    
    matResultsCluster = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenAll),1);...
        ones(length(cluster2withinAll),1)*2;ones(length(cluster2betweenAll),1)*2;...
        ones(length(cluster3withinAll),1)*3;ones(length(cluster3betweenAll),1)*3;...
        ones(length(cluster4withinAll),1)*4; ones(length(cluster4betweenAll),1)*4];
    
    matResultsType = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenAll),1);...
        ones(length(cluster2withinAll),1);ones(length(cluster2betweenAll),1);...
        ones(length(cluster3withinAll),1);ones(length(cluster3betweenAll),1);...
        ones(length(cluster4withinAll),1); ones(length(cluster4betweenAll),1)];
    
    matResultsBW = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenAll),1)*2;...
        ones(length(cluster2withinAll),1);ones(length(cluster2betweenAll),1)*2;...
        ones(length(cluster3withinAll),1);ones(length(cluster3betweenAll),1)*2;...
        ones(length(cluster4withinAll),1); ones(length(cluster4betweenAll),1)*2];
    
    matAll = [matResultsV, matResultsCluster, matResultsType, matResultsBW];
    
    f = figure; hold on;
    title('Linear regression slop (b1)');
    groupsBox1 = {};
    groupsBox1(end+1:end+length(indexLocation)) = {'Cluster1'};
    groupsBox1(end+1:end+length(indexLocation)) = {'Cluster2'};
    groupsBox1(end+1:end+length(indexLocation)) = {'Cluster3'};
    groupsBox1(end+1:end+length(indexLocation)) = {'Cluster4'};
    meanAll = [];
    bC = boxchart(categorical(groupsBox1), [tableResults.b1_cluster1(indexLocation); ...
        tableResults.b1_cluster2(indexLocation);...
        tableResults.b1_cluster3(indexLocation);...
        tableResults.b1_cluster4(indexLocation)]);
    meanAll(1) = mean(tableResults.b1_cluster1(indexLocation));
    meanAll(2) = mean(tableResults.b1_cluster2(indexLocation));
    meanAll(3) = mean(tableResults.b1_cluster3(indexLocation));   
    meanAll(4) = mean(tableResults.b1_cluster4(indexLocation));  
    plot(meanAll, '-*k');
    ylim([-1.2, 0.2]);
    mysave(f, fullfile(outputfolder, 'B1_slop'));   
    
    f = figure; hold on;
    title('Linear regression slop (b1) Hemi');
    groupsBox1 = {};
    groupsBox1(end+1:end+length(indexHemi)) = {'Cluster1'};
    groupsBox1(end+1:end+length(indexHemi)) = {'Cluster2'};
    groupsBox1(end+1:end+length(indexHemi)) = {'Cluster3'};
    groupsBox1(end+1:end+length(indexHemi)) = {'Cluster4'};
    meanAll = [];
    bC = boxchart(categorical(groupsBox1),...
        [min([tableResults.b1_sub1_cluster1(indexHemi),tableResults.b1_sub2_cluster1(indexHemi)],[],2); ...
        min([tableResults.b1_sub1_cluster2(indexHemi),tableResults.b1_sub2_cluster2(indexHemi)],[], 2);...
        min([tableResults.b1_sub1_cluster3(indexHemi),tableResults.b1_sub2_cluster3(indexHemi)],[], 2);...
        min([tableResults.b1_sub1_cluster4(indexHemi),tableResults.b1_sub2_cluster4(indexHemi)],[], 2)]);
    meanAll(1) = mean(min([tableResults.b1_sub1_cluster1(indexHemi),tableResults.b1_sub2_cluster1(indexHemi)],[],2));
    meanAll(2) = mean(min([tableResults.b1_sub1_cluster2(indexHemi),tableResults.b1_sub2_cluster2(indexHemi)],[],2));
    meanAll(3) = mean(min([tableResults.b1_sub1_cluster3(indexHemi),tableResults.b1_sub2_cluster3(indexHemi)],[],2));   
    meanAll(4) = mean(min([tableResults.b1_sub1_cluster4(indexHemi),tableResults.b1_sub2_cluster4(indexHemi)],[], 2));  
    plot(meanAll, '-*k');
    ylim([-1.2, 0.2]);
    mysave(f, fullfile(outputfolder, 'B1_slop_Hemi'));  
    
    statPlot([tableResults.b1_cluster1(indexLocation), ...
        tableResults.b1_cluster2(indexLocation),...
        tableResults.b1_cluster3(indexLocation),...
        tableResults.b1_cluster4(indexLocation)], 'slop(b1)', 'Linear regression slop (b1)', outputfolder, [0,0,128]./255);
    
    bslop = [tableResults.b1_cluster1(indexLocation), ...
        tableResults.b1_cluster2(indexLocation),...
        tableResults.b1_cluster3(indexLocation),...
        tableResults.b1_cluster4(indexLocation)];
end

function matAll = BetweenVsWithinSecDepth(indexLocation, tableResults, outputfolder, alphValue)
    cluster1betweenMain = [];
    cluster1betweenSec = [];
    cluster1withinAll = [];
    
    cluster2betweenMain = [];
    cluster2betweenSec = [];
    cluster2withinAll = [];
    
    cluster3betweenMain = [];
    cluster3betweenSec = [];
    cluster3withinAll = [];
    
    cluster4betweenMain = [];
    cluster4betweenSec = [];
    cluster4withinAll = [];


    for i = 1:size(indexLocation, 1)
        cluster1File = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\BetweenAndWithinSubTrees\AVSD*_Depth*.csv'];
        
        fileList = dir(cluster1File);
        
        if isempty(fileList)
            error('to add dir');
        end
        
        temp1 = [];
        temp1_1 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
                    
            if contains(fileList(j).name, '_Depth2_')
                if contains(fileList(j).name, 'BetweenMainDepth')
                    cluster1betweenMain(end+1) = nanmean(tableR.Var2);
                end
            elseif contains(fileList(j).name, '_BetweenSecondDepth') || contains(fileList(j).name, '_NotInDepth')
                temp1_1(end+1:end+length(tableR.Var2)) = tableR.Var2;
        
            else
                temp1(end+1:end+length(tableR.Var2)) = tableR.Var2;
        
            end
        end
        cluster1betweenSec(end+1) = nanmean(temp1_1);
            
        cluster1withinAll(end+1) = nanmean(temp1);
        
        cluster1withinAll(isnan(cluster1withinAll)) = [];
        cluster1betweenMain(isnan(cluster1betweenMain)) = [];
        cluster1betweenSec(isnan(cluster1betweenSec)) = [];
         
        cluster2File = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\BetweenAndWithinSubTrees\AVSD*_Depth*.csv'];
        
        fileList = dir(cluster2File);
        
        if isempty(fileList)
            error('to add dir');
        end
        
        temp2 = [];
        temp2_1 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
                    
            if contains(fileList(j).name, '_Depth2_')
                if contains(fileList(j).name, 'BetweenMainDepth')
                    cluster2betweenMain(end+1) = nanmean(tableR.Var2);
                end
            elseif contains(fileList(j).name, '_BetweenSecondDepth') || contains(fileList(j).name, '_NotInDepth')
                
                temp2_1(end+1:end+length(tableR.Var2)) = tableR.Var2;
            else
                temp2(end+1:end+length(tableR.Var2)) = tableR.Var2;
        
            end
        end
        
        
        cluster2betweenSec(end+1) = nanmean(temp2_1);
        cluster2withinAll(end+1) = nanmean(temp2);
        
        cluster2withinAll(isnan(cluster2withinAll)) = [];
        cluster2betweenMain(isnan(cluster2betweenMain)) = [];
        cluster2betweenSec(isnan(cluster2betweenSec)) = [];
         
        cluster3File = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\BetweenAndWithinSubTrees\AVSD*_Depth*.csv'];
        
        fileList = dir(cluster3File);
        
        if isempty(fileList)
            error('to add dir');
        end
        
        temp3 = [];
        temp3_1 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
                    
            if contains(fileList(j).name, '_Depth2_')
                if contains(fileList(j).name, 'BetweenMainDepth')
                    cluster3betweenMain(end+1) = nanmean(tableR.Var2);
                end
            elseif contains(fileList(j).name, '_BetweenSecondDepth') || contains(fileList(j).name, '_NotInDepth')
                temp3_1(end+1:end+length(tableR.Var2)) = tableR.Var2;
            else
                temp3(end+1:end+length(tableR.Var2)) = tableR.Var2;
        
            end
        end
        cluster3betweenSec(end+1) = nanmean(temp3_1);
            
        cluster3withinAll(end+1) = nanmean(temp3);
        
        cluster3withinAll(isnan(cluster3withinAll)) = [];
        cluster3betweenMain(isnan(cluster3betweenMain)) = [];
        cluster3betweenSec(isnan(cluster3betweenSec)) = [];
         
        cluster4File = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\BetweenAndWithinSubTrees\AVSD*_Depth*.csv'];
        
        fileList = dir(cluster4File);
        
        if isempty(fileList)
            error('to add dir');
        end
        
        temp4 = [];
        temp4_1 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
                    
            if contains(fileList(j).name, '_Depth2_')
                if contains(fileList(j).name, 'BetweenMainDepth')
                    cluster4betweenMain(end+1) = nanmean(tableR.Var2);
                end
            elseif contains(fileList(j).name, '_BetweenSecondDepth') || contains(fileList(j).name, '_NotInDepth')
                temp4_1(end+1:end+length(tableR.Var2)) = tableR.Var2;       
            else
                temp4(end+1:end+length(tableR.Var2)) = tableR.Var2;       
            end
        end
        cluster4betweenSec(end+1) = nanmean(temp4_1);
            
        cluster4withinAll(end+1) = nanmean(temp4);        
        
        cluster4withinAll(isnan(cluster4withinAll)) = [];
        cluster4betweenMain(isnan(cluster4betweenMain)) = [];
        cluster4betweenSec(isnan(cluster4betweenSec)) = [];     
    end
      
    f = figure; hold on;
    groupsBox1 = {};
    groupsBox1(end+1:end+length(cluster1withinAll)) = {'Within-Cluster1'};
    groupsBox1(end+1:end+length(cluster2withinAll)) = {'Within-Cluster2'};
    groupsBox1(end+1:end+length(cluster3withinAll)) = {'Within-Cluster3'};
    groupsBox1(end+1:end+length(cluster4withinAll)) = {'Within-Cluster4'};
    groupsBox2 = {};
    groupsBox2(end+1:end+length(cluster1betweenMain)) = {'Between-Main-Cluster1'};
    groupsBox2(end+1:end+length(cluster2betweenMain)) = {'Between-Main-Cluster2'};
    groupsBox2(end+1:end+length(cluster3betweenMain)) = {'Between-Main-Cluster3'};
    groupsBox2(end+1:end+length(cluster4betweenMain)) = {'Between-Main-Cluster4'};
    groupsBox3 = {};
    groupsBox3(end+1:end+length(cluster1betweenSec)) = {'Between-Sec-Cluster1'};
    groupsBox3(end+1:end+length(cluster2betweenSec)) = {'Between-Sec-Cluster2'};
    groupsBox3(end+1:end+length(cluster3betweenSec)) = {'Between-Sec-Cluster3'};
    groupsBox3(end+1:end+length(cluster4betweenSec)) = {'Between-Sec-Cluster4'};
    
    meanAll(1) = mean(cluster1withinAll);
    meanAll(2) = mean(cluster2withinAll);
    meanAll(3) = mean(cluster3withinAll);   
    meanAll(4) = mean(cluster4withinAll);  
    meanAll(5) = mean(cluster1betweenMain);
    meanAll(6) = mean(cluster2betweenMain);    
    meanAll(7) = mean(cluster3betweenMain);  
    meanAll(8) = mean(cluster4betweenMain);
    meanAll(9) = mean(cluster1betweenSec);
    meanAll(10) = mean(cluster2betweenSec);    
    meanAll(11) = mean(cluster3betweenSec);  
    meanAll(12) = mean(cluster4betweenSec);
    

    bC = boxchart(categorical(groupsBox1),[cluster1withinAll'; ...
        cluster2withinAll';...
        cluster3withinAll';...
        cluster4withinAll']);
    
    bC2 = boxchart(categorical(groupsBox2),[cluster1betweenMain';...
        cluster2betweenMain';...
        cluster3betweenMain';...
        cluster4betweenMain']);
    
    bC3 = boxchart(categorical(groupsBox3),[cluster1betweenSec';...
        cluster2betweenSec';...
        cluster3betweenSec';...
        cluster4betweenSec']);
    
    bC.BoxFaceColor = [0,0,0];
    bC.BoxFaceAlpha = 0.4;
    bC.MarkerColor = [0,0,0];
    bC2.BoxFaceColor = [0,0,255] ./ 255;
    bC2.MarkerColor = [0,0,255] ./ 255;
    bC2.BoxFaceAlpha = 0.4;
    bC3.BoxFaceColor = [0,255,0] ./ 255;
    bC3.MarkerColor = [0,255,0] ./ 255;
    bC3.BoxFaceAlpha = 0.4;
    plot(meanAll, '-*k');

    mysave(f, fullfile(outputfolder, 'PlotMeanAndSTDByDataWithinBetweenAllClusters_includeSub'));
    
    matResultsV = [cluster1withinAll';cluster1betweenMain';cluster1betweenSec';...
        cluster2withinAll';cluster2betweenMain';cluster2betweenSec';...
        cluster3withinAll';cluster3betweenMain';cluster3betweenSec';...
        cluster4withinAll';cluster4betweenMain';cluster4betweenSec'];
    
    matResultsCluster = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenMain),1); ones(length(cluster1betweenSec),1);...
        ones(length(cluster2withinAll),1)*2;ones(length(cluster2betweenMain),1)*2;ones(length(cluster2betweenSec),1)*2;...
        ones(length(cluster3withinAll),1)*3;ones(length(cluster3betweenMain),1)*3;ones(length(cluster3betweenSec),1)*3;...
        ones(length(cluster4withinAll),1)*4; ones(length(cluster4betweenMain),1)*4; ones(length(cluster4betweenSec),1)*4];
    
    matResultsType = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenMain),1);ones(length(cluster1betweenSec),1);...
        ones(length(cluster2withinAll),1);ones(length(cluster2betweenMain),1);ones(length(cluster2betweenSec),1);...
        ones(length(cluster3withinAll),1);ones(length(cluster3betweenMain),1);ones(length(cluster3betweenSec),1);...
        ones(length(cluster4withinAll),1); ones(length(cluster4betweenMain),1);ones(length(cluster4betweenSec),1)];
    
    matResultsBW = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenMain),1)*2;ones(length(cluster1betweenSec),1)*3;...
        ones(length(cluster2withinAll),1);ones(length(cluster2betweenMain),1)*2;ones(length(cluster2betweenSec),1)*3;...
        ones(length(cluster3withinAll),1);ones(length(cluster3betweenMain),1)*2;ones(length(cluster3betweenSec),1)*3;...
        ones(length(cluster4withinAll),1);ones(length(cluster4betweenMain),1)*2;ones(length(cluster4betweenSec),1)*3];
    
    matAll = [matResultsV, matResultsCluster, matResultsType, matResultsBW];
    save(fullfile(outputfolder, 'Results_MeanAndSTDByDataWithinBetweenAllClusters_includeSub.mat'), 'matAll');
end

function scatterPlotDataC(figT, outputfolder, temp2, temp2_x, temp2_1, temp2_1_x, temp2_2, temp2_2_x , temp2_3, temp2_3_x)
        f = figure; hold on;
        title(figT);
      
        colorD = [0.82,0.22,0.31];
        scatter(temp2_3_x ./ max(temp2_3_x), temp2_3, 'MarkerEdgeColor', colorD, 'MarkerFaceColor', colorD,'MarkerFaceAlpha' ,0.3, 'SizeData', 20);
        
        colorD = [0,130,200]./255;
        scatter(temp2_2_x ./ max(temp2_3_x), temp2_2, 'MarkerEdgeColor', colorD,'MarkerFaceColor', colorD,'MarkerFaceAlpha' ,0.3, 'SizeData', 20);
        
        colorD = [145,30,180] ./ 255;
        scatter(temp2_1_x ./ max(temp2_3_x), temp2_1,'MarkerEdgeColor', colorD, 'MarkerFaceColor', colorD,'MarkerFaceAlpha' ,0.3, 'SizeData', 20);
        
        colorD = [70, 240, 240] ./ 255;
        scatter(temp2_x ./ max(temp2_3_x), temp2,'MarkerEdgeColor', colorD, 'MarkerFaceColor', colorD,'MarkerFaceAlpha' ,0.3, 'SizeData', 20);
        
%         colorD = [0.82,0.22,0.31];
        colorD = [0, 0, 0] ./ 255;
        errorbar(mean(temp2_3_x ./ max(temp2_3_x)), mean(temp2_3), std(temp2_3), '-o','color', [0,0,0],'LineWidth', 1.2, 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', colorD, 'MarkerSize', 6);
        
        colorD = [0, 0, 0] ./ 255;
%         colorD = [145,30,180]./255;
        errorbar(mean(temp2_2_x ./ max(temp2_3_x)), mean(temp2_2), std(temp2_2), '-o','color', [0,0,0],'LineWidth', 1.2, 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', colorD, 'MarkerSize', 6);
       
        colorD = [0, 0, 0] ./ 255;
%         colorD = [245,130,48] ./ 255;
        errorbar(mean(temp2_1_x ./ max(temp2_3_x)), mean(temp2_1), std(temp2_1), '-o','color', [0,0,0],'LineWidth', 1.2, 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', colorD, 'MarkerSize', 6);
        
        colorD = [0, 0, 0] ./ 255;
%         colorD = [60, 180, 75] ./ 255;
        errorbar(mean(temp2_x ./ max(temp2_3_x)), mean(temp2), std(temp2), '-o','color', [0,0,0],'LineWidth', 1.2, 'MarkerEdgeColor', [0,0,0], 'MarkerFaceColor', colorD, 'MarkerSize', 6);
        
        ylim([-0.5,1]);
        mysave(f, fullfile(outputfolder, [figT, '_scatter_Plot']));
        
end


function matAll = BetweenVsWithinHemiDepth(indexLocation, tableResults, outputfolder, alphValue)
    cluster1betweenMain = [];
    cluster1betweenSec = [];
    cluster1withinAll = [];
    
    cluster2betweenMain = [];
    cluster2betweenSec = [];
    cluster2withinAll = [];
    
    cluster3betweenMain = [];
    cluster3betweenSec = [];
    cluster3withinAll = [];
    
    cluster4betweenMain = [];
    cluster4betweenSec = [];
    cluster4withinAll = [];

    cluster1betweenHemi = [];
    cluster2betweenHemi = [];
    cluster3betweenHemi = [];
    cluster4betweenHemi = [];

    if size(indexLocation,1) == 0
        matAll = nan(0,4);
        return;
    end
        
    for i = 1:size(indexLocation, 1)
        cluster1File = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\BetweenAndWithinSubTrees\AVSD*_Depth*.csv'];
        
        fileList = dir(cluster1File);
        
        if isempty(fileList)
            error('to add dir');
        end
        
        temp1 = [];
        temp1_1 = [];
        temp1_2 = [];
        temp1_3 = [];
        
        temp1_x = [];
        temp1_1_x = [];
        temp1_2_x = [];
        temp1_3_x= [];
        
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
                    
            if contains(fileList(j).name, '_Depth2_')
                if contains(fileList(j).name, 'BetweenMainDepth')
                    temp1_3(end+1:end+length(tableR.Var2)) = tableR.Var2; 
                    temp1_3_x(end+1:end+length(tableR.Var1)) = tableR.Var1;   
                end
            elseif contains(fileList(j).name, '_BetweenSecondDepth') || contains(fileList(j).name, '_NotInDepth')
                temp1_1(end+1:end+length(tableR.Var2)) = tableR.Var2;      
                temp1_1_x(end+1:end+length(tableR.Var1)) = tableR.Var1;  
            elseif contains(fileList(j).name, '_BetweenHemiDepth')
                temp1_2(end+1:end+length(tableR.Var2)) = tableR.Var2; 
                temp1_2_x(end+1:end+length(tableR.Var1)) = tableR.Var1;       
            else
                temp1(end+1:end+length(tableR.Var2)) = tableR.Var2; 
                temp1_x(end+1:end+length(tableR.Var1)) = tableR.Var1;      
            end
        end
        
        scatterPlotDataC('cluster 1', outputfolder, temp1, temp1_x, temp1_1, temp1_1_x, temp1_2, temp1_2_x , temp1_3, temp1_3_x)
        
        cluster1betweenMain(end+1) = nanmean(temp1_3);
        cluster1betweenSec(end+1) = nanmean(temp1_1);
        cluster1betweenHemi(end+1) = nanmean(temp1_2);
            
        cluster1withinAll(end+1) = nanmean(temp1);
        
        cluster1withinAll(isnan(cluster1withinAll)) = [];
        cluster1betweenMain(isnan(cluster1betweenMain)) = [];
        cluster1betweenSec(isnan(cluster1betweenSec)) = [];
        cluster1betweenHemi(isnan(cluster1betweenHemi)) = [];
         
        cluster2File = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\BetweenAndWithinSubTrees\AVSD*_Depth*.csv'];
        
        fileList = dir(cluster2File);
        
        if isempty(fileList)
            error('to add dir');
        end
        
        temp2 = [];
        temp2_1 = [];
        temp2_2 = [];
        temp2_3 = [];
        
        temp2_x = [];
        temp2_1_x = [];
        temp2_2_x = [];
        temp2_3_x = [];
        
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
                    
            if contains(fileList(j).name, '_Depth2_')
                if contains(fileList(j).name, 'BetweenMainDepth')
                    temp2_3(end+1:end+length(tableR.Var2)) = tableR.Var2;
                    temp2_3_x(end+1:end+length(tableR.Var1)) = tableR.Var1;
                end
            elseif contains(fileList(j).name, '_BetweenSecondDepth') || contains(fileList(j).name, '_NotInDepth')
                
                temp2_1(end+1:end+length(tableR.Var2)) = tableR.Var2;
                temp2_1_x(end+1:end+length(tableR.Var1)) = tableR.Var1;
            elseif contains(fileList(j).name, '_BetweenHemiDepth')
                temp2_2(end+1:end+length(tableR.Var2)) = tableR.Var2; 
                temp2_2_x(end+1:end+length(tableR.Var1)) = tableR.Var1;       
            else
                temp2(end+1:end+length(tableR.Var2)) = tableR.Var2;
                temp2_x(end+1:end+length(tableR.Var1)) = tableR.Var1;
        
            end
        end
        
        scatterPlotDataC('cluster 2', outputfolder, temp2, temp2_x, temp2_1, temp2_1_x, temp2_2, temp2_2_x , temp2_3, temp2_3_x)
        
        cluster2betweenMain(end+1) = nanmean(temp2_3);
        cluster2betweenHemi(end+1) = nanmean(temp2_2);
        cluster2betweenSec(end+1) = nanmean(temp2_1);
        cluster2withinAll(end+1) = nanmean(temp2);
        
        cluster2withinAll(isnan(cluster2withinAll)) = [];
        cluster2betweenMain(isnan(cluster2betweenMain)) = [];
        cluster2betweenSec(isnan(cluster2betweenSec)) = [];
        cluster2betweenHemi(isnan(cluster2betweenHemi)) = [];
        
        cluster3File = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\BetweenAndWithinSubTrees\AVSD*_Depth*.csv'];
        
        fileList = dir(cluster3File);
        
        if isempty(fileList)
            error('to add dir');
        end
        
        temp3 = [];
        temp3_1 = [];
        temp3_2 = [];
        temp3_3 = [];
        
        temp3_x = [];
        temp3_1_x = [];
        temp3_2_x = [];
        temp3_3_x = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
                    
            if contains(fileList(j).name, '_Depth2_')
                if contains(fileList(j).name, 'BetweenMainDepth')
                    temp3_3(end+1:end+length(tableR.Var2)) = tableR.Var2;
                    temp3_3_x(end+1:end+length(tableR.Var1)) = tableR.Var1;
                end
            elseif contains(fileList(j).name, '_BetweenSecondDepth') || contains(fileList(j).name, '_NotInDepth')
                temp3_1(end+1:end+length(tableR.Var2)) = tableR.Var2;
                temp3_1_x(end+1:end+length(tableR.Var1)) = tableR.Var1;
            elseif contains(fileList(j).name, '_BetweenHemiDepth')
                temp3_2(end+1:end+length(tableR.Var2)) = tableR.Var2;
                temp3_2_x(end+1:end+length(tableR.Var1)) = tableR.Var1;        
            else
                temp3(end+1:end+length(tableR.Var2)) = tableR.Var2;
                temp3_x(end+1:end+length(tableR.Var1)) = tableR.Var1;
        
            end
        end
        
        scatterPlotDataC('cluster 3', outputfolder, temp3, temp3_x, temp3_1, temp3_1_x, temp3_2, temp3_2_x , temp3_3, temp3_3_x)
        
        cluster3betweenMain(end+1) = nanmean(temp3_3);
        cluster3betweenSec(end+1) = nanmean(temp3_1);
        cluster3betweenHemi(end+1) = nanmean(temp3_2);
            
        cluster3withinAll(end+1) = nanmean(temp3);
        
        cluster3withinAll(isnan(cluster3withinAll)) = [];
        cluster3betweenMain(isnan(cluster3betweenMain)) = [];
        cluster3betweenSec(isnan(cluster3betweenSec)) = [];
        cluster3betweenHemi(isnan(cluster3betweenHemi)) = [];
         
        cluster4File = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\BetweenAndWithinSubTrees\AVSD*_Depth*.csv'];
        
        fileList = dir(cluster4File);
        
        if isempty(fileList)
            error('to add dir');
        end
        
        temp4 = [];
        temp4_1 = [];
        temp4_2 = [];
        temp4_3 = [];
        
        temp4_x = [];
        temp4_1_x = [];
        temp4_2_x = [];
        temp4_3_x = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
                    
            if contains(fileList(j).name, '_Depth2_')
                if contains(fileList(j).name, 'BetweenMainDepth')
                     temp4_3(end+1:end+length(tableR.Var2)) = tableR.Var2; 
                     temp4_3_x(end+1:end+length(tableR.Var1)) = tableR.Var1;   
                end
            elseif contains(fileList(j).name, '_BetweenSecondDepth') || contains(fileList(j).name, '_NotInDepth')
                temp4_1(end+1:end+length(tableR.Var2)) = tableR.Var2;  
                temp4_1_x(end+1:end+length(tableR.Var1)) = tableR.Var1;     
            elseif contains(fileList(j).name, '_BetweenHemiDepth')
                temp4_2(end+1:end+length(tableR.Var2)) = tableR.Var2;
                temp4_2_x(end+1:end+length(tableR.Var1)) = tableR.Var1;        
            else
                temp4(end+1:end+length(tableR.Var2)) = tableR.Var2;
                temp4_x(end+1:end+length(tableR.Var1)) = tableR.Var1;       
            end
        end
        
        scatterPlotDataC('cluster 4', outputfolder, temp4, temp4_x, temp4_1, temp4_1_x, temp4_2, temp4_2_x , temp4_3, temp4_3_x)
        
        cluster4betweenMain(end+1) = nanmean(temp4_3);
        cluster4betweenSec(end+1) = nanmean(temp4_1);
        cluster4betweenHemi(end+1) = nanmean(temp4_2);
            
        cluster4withinAll(end+1) = nanmean(temp4);        
        
        cluster4withinAll(isnan(cluster4withinAll)) = [];
        cluster4betweenMain(isnan(cluster4betweenMain)) = [];
        cluster4betweenSec(isnan(cluster4betweenSec)) = []; 
        cluster4betweenHemi(isnan(cluster4betweenSec)) = [];     
    end
      
    f = figure; hold on;
    groupsBox1 = {};
    groupsBox1(end+1:end+length(cluster1withinAll)) = {'Within-Cluster1'};
    groupsBox1(end+1:end+length(cluster2withinAll)) = {'Within-Cluster2'};
    groupsBox1(end+1:end+length(cluster3withinAll)) = {'Within-Cluster3'};
    groupsBox1(end+1:end+length(cluster4withinAll)) = {'Within-Cluster4'};
    groupsBox2 = {};
    groupsBox2(end+1:end+length(cluster1betweenMain)) = {'Between-Main-Cluster1'};
    groupsBox2(end+1:end+length(cluster2betweenMain)) = {'Between-Main-Cluster2'};
    groupsBox2(end+1:end+length(cluster3betweenMain)) = {'Between-Main-Cluster3'};
    groupsBox2(end+1:end+length(cluster4betweenMain)) = {'Between-Main-Cluster4'};
    groupsBox3 = {};
    groupsBox3(end+1:end+length(cluster1betweenSec)) = {'Between-Sec-Cluster1'};
    groupsBox3(end+1:end+length(cluster2betweenSec)) = {'Between-Sec-Cluster2'};
    groupsBox3(end+1:end+length(cluster3betweenSec)) = {'Between-Sec-Cluster3'};
    groupsBox3(end+1:end+length(cluster4betweenSec)) = {'Between-Sec-Cluster4'};
    
    meanAll(1) = mean(cluster1withinAll);
    meanAll(2) = mean(cluster2withinAll);
    meanAll(3) = mean(cluster3withinAll);   
    meanAll(4) = mean(cluster4withinAll);  
    meanAll(5) = mean(cluster1betweenMain);
    meanAll(6) = mean(cluster2betweenMain);    
    meanAll(7) = mean(cluster3betweenMain);  
    meanAll(8) = mean(cluster4betweenMain);
    meanAll(9) = mean(cluster1betweenSec);
    meanAll(10) = mean(cluster2betweenSec);    
    meanAll(11) = mean(cluster3betweenSec);  
    meanAll(12) = mean(cluster4betweenSec);
    

    bC = boxchart(categorical(groupsBox1),[cluster1withinAll'; ...
        cluster2withinAll';...
        cluster3withinAll';...
        cluster4withinAll']);
    
    bC2 = boxchart(categorical(groupsBox2),[cluster1betweenMain';...
        cluster2betweenMain';...
        cluster3betweenMain';...
        cluster4betweenMain']);
    
    bC3 = boxchart(categorical(groupsBox3),[cluster1betweenSec';...
        cluster2betweenSec';...
        cluster3betweenSec';...
        cluster4betweenSec']);
    
    bC.BoxFaceColor = [0,0,0];
    bC.BoxFaceAlpha = 0.4;
    bC.MarkerColor = [0,0,0];
    bC2.BoxFaceColor = [0,0,255] ./ 255;
    bC2.MarkerColor = [0,0,255] ./ 255;
    bC2.BoxFaceAlpha = 0.4;
    bC3.BoxFaceColor = [0,255,0] ./ 255;
    bC3.MarkerColor = [0,255,0] ./ 255;
    bC3.BoxFaceAlpha = 0.4;
    plot(meanAll, '-*k');

    mysave(f, fullfile(outputfolder, 'PlotMeanAndSTDByDataWithinBetweenAllClusters_includeSub'));
    
    matResultsV = [cluster1withinAll';cluster1betweenMain';cluster1betweenSec';...
        cluster2withinAll';cluster2betweenMain';cluster2betweenSec';...
        cluster3withinAll';cluster3betweenMain';cluster3betweenSec';...
        cluster4withinAll';cluster4betweenMain';cluster4betweenSec'];
    
    matResultsCluster = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenMain),1); ones(length(cluster1betweenSec),1);...
        ones(length(cluster2withinAll),1)*2;ones(length(cluster2betweenMain),1)*2;ones(length(cluster2betweenSec),1)*2;...
        ones(length(cluster3withinAll),1)*3;ones(length(cluster3betweenMain),1)*3;ones(length(cluster3betweenSec),1)*3;...
        ones(length(cluster4withinAll),1)*4; ones(length(cluster4betweenMain),1)*4; ones(length(cluster4betweenSec),1)*4];
    
    matResultsType = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenMain),1);ones(length(cluster1betweenSec),1);...
        ones(length(cluster2withinAll),1);ones(length(cluster2betweenMain),1);ones(length(cluster2betweenSec),1);...
        ones(length(cluster3withinAll),1);ones(length(cluster3betweenMain),1);ones(length(cluster3betweenSec),1);...
        ones(length(cluster4withinAll),1); ones(length(cluster4betweenMain),1);ones(length(cluster4betweenSec),1)];
    
    matResultsBW = [ones(length(cluster1withinAll),1); ones(length(cluster1betweenMain),1)*2;ones(length(cluster1betweenSec),1)*3;...
        ones(length(cluster2withinAll),1);ones(length(cluster2betweenMain),1)*2;ones(length(cluster2betweenSec),1)*3;...
        ones(length(cluster3withinAll),1);ones(length(cluster3betweenMain),1)*2;ones(length(cluster3betweenSec),1)*3;...
        ones(length(cluster4withinAll),1);ones(length(cluster4betweenMain),1)*2;ones(length(cluster4betweenSec),1)*3];
    
    matAll = [matResultsV, matResultsCluster, matResultsType, matResultsBW];
    save(fullfile(outputfolder, 'Results_MeanAndSTDByDataWithinBetweenAllClusters_includeSub.mat'), 'matAll');
    
%     cluster1withinAll(cluster1withinAll < 0 ) = 0;
%     cluster2withinAll(cluster2withinAll < 0 ) = 0;
%     cluster3withinAll(cluster3withinAll < 0 ) = 0;
%     cluster4withinAll(cluster4withinAll < 0 ) = 0;
%     cluster1betweenHemi(cluster1betweenHemi < 0 ) = 0;
%     cluster2betweenHemi(cluster2betweenHemi < 0 ) = 0;
%     cluster3betweenHemi(cluster3betweenHemi < 0 ) = 0;
%     cluster4betweenHemi(cluster4betweenHemi < 0 ) = 0;
%     cluster1betweenSec(cluster1betweenSec < 0 ) = 0;
%     cluster2betweenSec(cluster2betweenSec < 0 ) = 0;
%     cluster3betweenSec(cluster3betweenSec < 0 ) = 0;
%     cluster4betweenSec(cluster4betweenSec < 0 ) = 0;
%     cluster1betweenMain(cluster1betweenMain < 0 ) = 0;
%     cluster2betweenMain(cluster2betweenMain < 0 ) = 0;
%     cluster3betweenMain(cluster3betweenMain < 0 ) = 0;
%     cluster4betweenMain(cluster4betweenMain < 0 ) = 0;
%     
    listN = {'Within quaternary bifurcations','Between quaternary bifurcations', 'Between tertiary bifurcations', 'Between secondary bifurcations'};
    f = figure; hold on;
    title('Cluster 1');
    bar(categorical(listN, listN), [mean(cluster1withinAll),mean(cluster1betweenHemi), mean(cluster1betweenSec), mean(cluster1betweenMain)]);
    ylim([-0.3,1]);
    mysave(f, fullfile(outputfolder, 'Cluster1_cell_3compare'));
    f = figure; hold on;
    title('Cluster 2');
    bar(categorical(listN, listN), [mean(cluster2withinAll),mean(cluster2betweenHemi), mean(cluster2betweenSec), mean(cluster2betweenMain)]);
    ylim([-0.3,1]);
    mysave(f, fullfile(outputfolder, 'Cluster2_cell_3compare'));
    f = figure; hold on;
    title('Cluster 3');
    bar(categorical(listN, listN), [mean(cluster3withinAll),mean(cluster3betweenHemi), mean(cluster3betweenSec), mean(cluster3betweenMain)]);
    ylim([-0.3,1]);
    mysave(f, fullfile(outputfolder, 'Cluster3_cell_3compare'));
    f = figure; hold on;
    title('Cluster 4');
    bar(categorical(listN, listN), [mean(cluster4withinAll),mean(cluster4betweenHemi), mean(cluster4betweenSec), mean(cluster4betweenMain)]);
    ylim([-0.3,1]);
    mysave(f, fullfile(outputfolder, 'Cluster4_cell_3compare'));
end

function statPlot(arrayT, ylabelName, TitleV, outputfolder, colorType)
    f = figure; hold on;
    b = boxchart(arrayT);
    ylim([0,1]);
    xticklabels({'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4'});
    b.BoxFaceColor = colorType;
    b.BoxFaceAlpha = 0.4;
    b.MarkerColor = colorType;
    ylabel(ylabelName);
    plot(mean(arrayT), '-*k');
    mysave(f, fullfile(outputfolder, TitleV));
    mResult = mean(arrayT);
    stdR = std(arrayT);
    textAnova = '';

    textAnova = strcat(textAnova, sprintf('cluster1 mean: %.4f, std: %.4f \\n ', mResult(1), stdR(1)));
    textAnova = strcat(textAnova, sprintf('cluster2 mean: %.4f, std: %.4f \\n ', mResult(2), stdR(2)));
    textAnova = strcat(textAnova, sprintf('cluster3 mean: %.4f, std: %.4f \\n ', mResult(3), stdR(3)));
    textAnova = strcat(textAnova, sprintf('cluster4 mean: %.4f, std: %.4f \\n ', mResult(4), stdR(4)));             

    fid=fopen(fullfile(outputfolder, ['summaryMandstd_', TitleV,'.txt']),'w');
    fprintf(fid, textAnova);
    fclose(fid);
    
    textAnova = '';
    [p,~,statsM] = anova1(arrayT, {'Cluster1', 'Cluster2', 'Cluster3', 'Cluster4'}, 'off');
    if p < 0.05
        [c,~,~,groupnames] = multcompare(statsM, 'Display', 'off');
        for j = 1:size(c, 1)
            textAnova = strcat(textAnova, sprintf('%s vs %s - pValue : %.4f, mean Diff : %.4f, CI_L : %.4f, CI_H : %.4f \\n ', ...
                groupnames{c(j, 1)}, groupnames{c(j, 2)}, c(j, 6), c(j, 4), c(j, 3), c(j, 5)));
        end
    end

    fid=fopen(fullfile(outputfolder, ['statistic_', TitleV,'.txt']),'w');
    fprintf(fid, textAnova);
    fclose(fid);
end