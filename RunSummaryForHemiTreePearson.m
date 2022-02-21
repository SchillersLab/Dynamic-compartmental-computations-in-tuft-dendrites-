function RunSummaryForHemiTreePearson()
    mantelRunnerLocation = '\\jackie-analysis\e\Shay\RunnersLocationSummary.xlsx';
    sheetName = 'RunOnlyTuft';
    outputfolder(1) = {'\\jackie-analysis\e\Shay\StatisticSummary\HemiTreePostSummary\Big_Fix5\'};
    mkdir(outputfolder{1});
    classification(1) = 1;
    
    outputfolder(2) = {'\\jackie-analysis\e\Shay\StatisticSummary\HemiTreePostSummary\Small_Fix5\'};
    mkdir(outputfolder{2});
    classification(2) = 2;
    
    colorsPerType = zeros([2,3]);
    colorsPerType(2, :) = [0,0,1];
    
    pcaMeanAll = {};
    pcaChance = {};
    pcaSTDAll = {};
    matAll = [];
    tableResults = readtable(mantelRunnerLocation,'Sheet',sheetName);
    indexLocation1 = find(tableResults.Classification == classification(1) & tableResults.R2HemiTreeInclude == 1);
    indexLocation2 = find(tableResults.Classification == classification(2) & tableResults.R2HemiTreeInclude == 1);
       
    for k = 1:2
        indexLocation = find(tableResults.Classification == classification(k) & tableResults.R2HemiTreeInclude == 1);
        indexLocationMain = find(tableResults.Classification == classification(k) & tableResults.includeMainR2 == 1);
       
        
        matRTemp2 = BetweenVsWithinSecDepth(indexLocationMain, tableResults, outputfolder{k}, 0.05);
        matRTemp2(:, 3) = k;
        
        matAll2(end+1:end + size(matRTemp2, 1),1:4) = matRTemp2;
        
        
        matRTemp = BetweenVsWithin(indexLocationMain, tableResults, outputfolder{k}, 0.05);
        matRTemp(:, 3) = k;
        
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
    title('PCA accuracy Mean')
    groupsBox1 = {};
    groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster1'};
    groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster2'};
    groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster3'};
    groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster4'};
    groupsBox2 = {};
    groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster1'};
    groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster2'};
    groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster3'};
    groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster4'};

    meanAll(1:4) = mean(pcaMeanAll{1});
    meanAll(5:8) = mean(pcaMeanAll{2});
    
    chaAll(1:4) = ones(1,4)*(pcaChance{1});
    chaAll(5:8) = ones(1,4)*(pcaChance{2});
    


    bC = boxchart(categorical(groupsBox1),[pcaMeanAll{1}(:,1); pcaMeanAll{1}(:,2); pcaMeanAll{1}(:,3); pcaMeanAll{1}(:,4)]);
    
    bC2 = boxchart(categorical(groupsBox2),[pcaMeanAll{2}(:,1); pcaMeanAll{2}(:,2); pcaMeanAll{2}(:,3); pcaMeanAll{2}(:,4)]);
    
  
    bC.BoxFaceColor = [0,0,0];
    bC.BoxFaceAlpha = 0.4;
    bC.MarkerColor = [0,0,0];
    bC2.BoxFaceColor = [0,0,255] ./ 255;
    bC2.MarkerColor = [0,0,255] ./ 255;
    bC2.BoxFaceAlpha = 0.4;
    
    plot(meanAll, '-*k');
    plot(chaAll, '--k')
    
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
    groupsBox1 = {};
    groupsBox1(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation1))) = {'Big-Cluster1'};
    groupsBox1(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation1))) = {'Big-Cluster2'};
    groupsBox1(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation1))) = {'Big-Cluster3'};
    groupsBox1(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation1))) = {'Big-Cluster4'};
    groupsBox2 = {};
    groupsBox2(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation2))) = {'Small-Cluster1'};
    groupsBox2(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation2))) = {'Small-Cluster2'};
    groupsBox2(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation2))) = {'Small-Cluster3'};
    groupsBox2(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation2))) = {'Small-Cluster4'};
    
    meanAll(1) = mean(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation1));
    meanAll(2) = mean(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation1));
    meanAll(3) = mean(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation1));   
    meanAll(4) = mean(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation1));  
    meanAll(5) = mean(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation2));
    meanAll(6) = mean(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation2));    
    meanAll(7) = mean(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation2));  
    meanAll(8) = mean(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation2));
    

    bC = boxchart(categorical(groupsBox1),[tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation1); ...
        tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation1);...
        tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation1);...
        tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation1)]);
    
    bC2 = boxchart(categorical(groupsBox2),[tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation2);...
        tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation2);...
        tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation2);...
        tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation2)]);
    
  
    bC.BoxFaceColor = [0,0,0];
    bC.BoxFaceAlpha = 0.4;
    bC.MarkerColor = [0,0,0];
    bC2.BoxFaceColor = [0,0,255] ./ 255;
    bC2.MarkerColor = [0,0,255] ./ 255;
    bC2.BoxFaceAlpha = 0.4;
    
    plot(meanAll, '-*k');
    
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
    calcAnovaForPearsonCompare2(matAll2, outputfolder{1});
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
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Big = %f \\n ',i, mean(y1)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Big = %f \\n ',i, mean(y2)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Small = %f \\n ',i, mean(y3)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Small = %f \\n ',i, mean(y4)));
        
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
        meanAll1 = [mean(y1), mean(y2);mean(y5)];
        meanAll2 = [mean(y3), mean(y4);mean(y6)];
        
        c1 = cell(1, length(y1));
        c1(:) = {'Big-Within'};
        c2 = cell(1, length(y2));
        c2(:) = {'Big-BetweenMain'};
        c5 = cell(1, length(y5));
        c5(:) = {'Big-BetweenSec'};
        c3 = cell(1, length(y3));
        c3(:) = {'Small-Within'};
        c4 = cell(1, length(y4));
        c4(:) = {'Small-BetweenMain'};
        c6 = cell(1, length(y6));
        c6(:) = {'Small-BetweenSec'};
        
        bc1 = boxchart(categorical(c1),y1);
        bc2 = boxchart(categorical(c2),y2);
        bc5 = boxchart(categorical(c5),y5);
        bc3 = boxchart(categorical(c3),y3);
        bc4 = boxchart(categorical(c4),y4);
        bc6 = boxchart(categorical(c6),y6);
        
        plot(categorical({'Big-Within', 'Big-Between-Main', 'Big-Between-Sec'}), meanAll1, '-*k');
        plot(categorical({'Small-Within', 'Small-Between', 'Small-Sec'}), meanAll2, '-*k');
        
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


function matAll = BetweenVsWithin(indexLocation, tableResults, outputfolder, alphValue)
    cluster1betweenAll = [];
    cluster1withinAll = [];
    cluster2betweenAll = [];
    cluster2withinAll = [];
    cluster3betweenAll = [];
    cluster3withinAll = [];
    cluster4betweenAll = [];
    cluster4withinAll = [];


    for i = 1:size(indexLocation, 1)
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*bp*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\AVSD*_Depth2_*bp*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster1betweenAll(end+1) = nanmean(tableBetween.Var2);
        xVals = [];
        xVals(end+1:end+length(tableBetween.Var1)) = tableBetween.Var1;
        yVals = [];
        yVals(end+1:end+length(tableBetween.Var2)) = tableBetween.Var2;
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster1withinAll(end+1) = nanmean(tablewithin.Var2);
            
            xVals(end+1:end+length(tablewithin.Var1)) = tablewithin.Var1;
            yVals(end+1:end+length(tablewithin.Var2)) = tablewithin.Var2;
        end
        
        cluster1withinAll(isnan(cluster1withinAll)) = [];
        cluster1betweenAll(isnan(cluster1betweenAll)) = [];
         
        mdAll1 = fitglm(xVals, yVals);
        
        
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*bp*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\AVSD*_Depth2_*bp*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster2betweenAll(end+1) = nanmean(tableBetween.Var2);
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster2withinAll(end+1) = nanmean(tablewithin.Var2);
        end
        
        cluster2withinAll(isnan(cluster2withinAll)) = [];
        cluster2betweenAll(isnan(cluster2betweenAll)) = [];
         
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*bp*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\AVSD*_Depth2_*bp*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster3betweenAll(end+1) = nanmean(tableBetween.Var2);
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster3withinAll(end+1) = nanmean(tablewithin.Var2);
        end
        
        cluster3withinAll(isnan(cluster3withinAll)) = [];
        cluster3betweenAll(isnan(cluster3betweenAll)) = [];
           
        
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*bp*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\AVSD*_Depth2_*bp*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster4betweenAll(end+1) = nanmean(tableBetween.Var2);
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster4withinAll(end+1) = nanmean(tablewithin.Var2);
        end
        
        cluster4withinAll(isnan(cluster4withinAll)) = [];
        cluster4betweenAll(isnan(cluster4betweenAll)) = [];
        
        
    end
    
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
        cluster1BetweenMainDepthFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1BetweenSecDepthFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\BetweenAndWithinSubTrees\AVSD*_Depth_*BetweenSecond.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\BetweenAndWithinSubTrees\AVSD*_Depth_*bp*.csv'];
        
        fileBetweenMain = dir(cluster1BetweenMainDepthFile);
        fileBetweenSec = dir(cluster1BetweenSecDepthFile);
        
        if isempty(fileBetweenMain)
            error('to add dir');
        end
        
        tableBetweenMain = readtable(fullfile(fileBetweenMain(1).folder,fileBetweenMain(1).name));
        tableBetweenSec = readtable(fullfile(fileBetweenSec(1).folder,fileBetweenSec(1).name));
        cluster1betweenMain(end+1) = nanmean(tableBetweenMain.Var2);
        cluster1betweenSec(end+1) = nanmean(tableBetweenSec.Var2);
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster1withinAll(end+1) = nanmean(tablewithin.Var2);
        end
        
        cluster1withinAll(isnan(cluster1withinAll)) = [];
        cluster1betweenMain(isnan(cluster1betweenMain)) = [];
        cluster1betweenSec(isnan(cluster1betweenSec)) = [];
         
        cluster2BetweenMainDepthFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster2BetweenSecDepthFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\BetweenAndWithinSubTrees\AVSD*_Depth_*BetweenSecond.csv'];
        cluster2WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\BetweenAndWithinSubTrees\AVSD*_Depth_*bp*.csv'];
        
        fileBetweenMain = dir(cluster2BetweenMainDepthFile);
        fileBetweenSec = dir(cluster2BetweenSecDepthFile);
        
        if isempty(fileBetweenMain)
            error('to add dir');
        end
        
        tableBetweenMain = readtable(fullfile(fileBetweenMain(1).folder,fileBetweenMain(1).name));
        tableBetweenSec = readtable(fullfile(fileBetweenSec(1).folder,fileBetweenSec(1).name));
        cluster2betweenMain(end+1) = nanmean(tableBetweenMain.Var2);
        cluster2betweenSec(end+1) = nanmean(tableBetweenSec.Var2);
        
        filewithin = dir(cluster2WithInFile);
        for j = 1:size(filewithin, 1)
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster2withinAll(end+1) = nanmean(tablewithin.Var2);
        end
        
        cluster2withinAll(isnan(cluster2withinAll)) = [];
        cluster2betweenMain(isnan(cluster2betweenMain)) = [];
        cluster2betweenSec(isnan(cluster2betweenSec)) = [];
         
        cluster3BetweenMainDepthFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster3BetweenSecDepthFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\BetweenAndWithinSubTrees\AVSD*_Depth_*BetweenSecond.csv'];
        cluster3WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\BetweenAndWithinSubTrees\AVSD*_Depth_*bp*.csv'];
        
        fileBetweenMain = dir(cluster3BetweenMainDepthFile);
        fileBetweenSec = dir(cluster3BetweenSecDepthFile);
        
        if isempty(fileBetweenMain)
            error('to add dir');
        end
        
        tableBetweenMain = readtable(fullfile(fileBetweenMain(1).folder,fileBetweenMain(1).name));
        tableBetweenSec = readtable(fullfile(fileBetweenSec(1).folder,fileBetweenSec(1).name));
        cluster3betweenMain(end+1) = nanmean(tableBetweenMain.Var2);
        cluster3betweenSec(end+1) = nanmean(tableBetweenSec.Var2);
        
        filewithin = dir(cluster3WithInFile);
        for j = 1:size(filewithin, 1)
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster3withinAll(end+1) = nanmean(tablewithin.Var2);
        end
        
        cluster3withinAll(isnan(cluster3withinAll)) = [];
        cluster3betweenMain(isnan(cluster3betweenMain)) = [];
        cluster3betweenSec(isnan(cluster3betweenSec)) = [];
         
        cluster4BetweenMainDepthFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster4BetweenSecDepthFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\BetweenAndWithinSubTrees\AVSD*_Depth_*BetweenSecond.csv'];
        cluster4WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\BetweenAndWithinSubTrees\AVSD*_Depth_*bp*.csv'];
        
        fileBetweenMain = dir(cluster4BetweenMainDepthFile);
        fileBetweenSec = dir(cluster4BetweenSecDepthFile);
        
        if isempty(fileBetweenMain)
            error('to add dir');
        end
        
        tableBetweenMain = readtable(fullfile(fileBetweenMain(1).folder,fileBetweenMain(1).name));
        tableBetweenSec = readtable(fullfile(fileBetweenSec(1).folder,fileBetweenSec(1).name));
        cluster4betweenMain(end+1) = nanmean(tableBetweenMain.Var2);
        cluster4betweenSec(end+1) = nanmean(tableBetweenSec.Var2);
        
        filewithin = dir(cluster4WithInFile);
        for j = 1:size(filewithin, 1)
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster4withinAll(end+1) = nanmean(tablewithin.Var2);
        end
        
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
    
    matResultsV = [cluster1withinAll';cluster1betweenAll';cluster2withinAll';cluster2betweenAll';...
        cluster3withinAll';cluster3betweenAll';cluster4withinAll';cluster4betweenAll'];
    
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
end

function BetweenVsWithinOLD(indexLocation, tableResults, outputfolder, alphValue)
    cluster1betweenAll = [];
    cluster1withinAll = [];
    cluster2betweenAll = [];
    cluster2withinAll = [];
    cluster3betweenAll = [];
    cluster3withinAll = [];
    cluster4betweenAll = [];
    cluster4withinAll = [];


    for i = 1:size(indexLocation, 1)
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\*\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\*\AVSD*_Depth2_*bp*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\AVSD*_Depth2_*bp*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster1betweenAll(end+1:end+length(tableBetween.Var2)) = tableBetween.Var2;
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster1withinAll(end+1:end+length(tablewithin.Var2)) = tablewithin.Var2;
        end
        
        cluster1withinAll(isnan(cluster1withinAll)) = [];
        cluster1betweenAll(isnan(cluster1betweenAll)) = [];
         
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\*\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\*\AVSD*_Depth2_*bp*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\AVSD*_Depth2_*bp*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster2betweenAll(end+1:end+length(tableBetween.Var2)) = tableBetween.Var2;
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster2withinAll(end+1:end+length(tablewithin.Var2)) = tablewithin.Var2;
        end
        
        cluster2withinAll(isnan(cluster2withinAll)) = [];
        cluster2betweenAll(isnan(cluster2betweenAll)) = [];
         
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\*\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\*\AVSD*_Depth2_*bp*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\AVSD*_Depth2_*bp*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster3betweenAll(end+1:end+length(tableBetween.Var2)) = tableBetween.Var2;
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster3withinAll(end+1:end+length(tablewithin.Var2)) = tablewithin.Var2;
        end
        
        cluster3withinAll(isnan(cluster3withinAll)) = [];
        cluster3betweenAll(isnan(cluster3betweenAll)) = [];
           
        
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\*\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\*\AVSD*_Depth2_*bp*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\AVSD*_Depth2_*bp*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster4betweenAll(end+1:end+length(tableBetween.Var2)) = tableBetween.Var2;
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster4withinAll(end+1:end+length(tablewithin.Var2)) = tablewithin.Var2;
        end
        
        cluster4withinAll(isnan(cluster4withinAll)) = [];
        cluster4betweenAll(isnan(cluster4betweenAll)) = [];
        
        
    end
    
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
    
%     meanAll = [nanmean(cluster1withinAll), nanmean(cluster1betweenAll),...
%         nanmean(cluster2withinAll), nanmean(cluster2betweenAll),...
%         nanmean(cluster3withinAll), nanmean(cluster3betweenAll),...
%         nanmean(cluster4withinAll), nanmean(cluster4betweenAll)];
%     stdAll = [std(cluster1withinAll), std(cluster1betweenAll),...
%         std(cluster2withinAll), std(cluster2betweenAll),...
%         std(cluster3withinAll), std(cluster3betweenAll),...
%         std(cluster4withinAll), std(cluster4betweenAll)];
    
    
    
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