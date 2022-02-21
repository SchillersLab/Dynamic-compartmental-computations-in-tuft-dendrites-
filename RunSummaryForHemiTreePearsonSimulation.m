function RunSummaryForHemiTreePearsonSimulation()
    mantelRunnerLocation = '\\jackie-analysis\e\Shay\RunnersLocationSummary.xlsx';
    sheetName = 'Simulation';
    outputfolder(1) = {'\\jackie-analysis\e\Shay\StatisticSummary\HemiTreePostSummary_Simulation\Big3\'};
    mkdir(outputfolder{1});
    classification(1) = 1;
    
    outputfolder(2) = {'\\jackie-analysis\e\Shay\StatisticSummary\HemiTreePostSummary_Simulation\SmallHZ3\'};
    mkdir(outputfolder{2});
    classification(2) = 2;
    
    outputfolder(3) = {'\\jackie-analysis\e\Shay\StatisticSummary\HemiTreePostSummary_Simulation\SmallHZTuft\'};
    mkdir(outputfolder{3});
    classification(3) = 3;
    
    colorsPerType = zeros([2,3]);
    colorsPerType(2, :) = [0,0,1];
    
    pcaMeanAll = {};
    pcaChance = {};
    pcaSTDAll = {};
    matAll = [];
    tableResults = readtable(mantelRunnerLocation,'Sheet',sheetName);
    indexLocation1 = find(tableResults.Classification == classification(1) & tableResults.R2HemiTreeInclude == 1);
    indexLocation2 = find(tableResults.Classification == classification(2) & tableResults.R2HemiTreeInclude == 1);
       
    for k = 1:3
        indexLocation = find(tableResults.Classification == classification(k) & tableResults.R2HemiTreeInclude == 1);
        indexLocationMain = find(tableResults.Classification == classification(k) & tableResults.includeMainR2 == 1);
       
        matRTemp = BetweenVsWithin(indexLocationMain, tableResults, outputfolder{k}, 0.05);
        matRTemp(:, 3) = k;
        
        matAll(end+1:end + size(matRTemp, 1),1:4) = matRTemp;
        
%         statPlot([tableResults.R2Cluster1_ALL(indexLocationMain),...
%             tableResults.R2Cluster2_ALL(indexLocationMain),...
%             tableResults.R2Cluster3_ALL(indexLocationMain),...
%             tableResults.R2Cluster4_ALL(indexLocationMain)], 'R2', 'R2Main', outputfolder{k}, colorsPerType(k, :))
%        
%         statPlot([tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation),...
%             tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation),...
%             tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation),...
%             tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation)], 'R2', 'R2maxHemiTree', outputfolder{k}, colorsPerType(k, :))
% 
%         statPlot([tableResults.R2Cluster1_HemiTree(indexLocation),...
%             tableResults.R2Cluster2_HemiTree(indexLocation),...
%             tableResults.R2Cluster3_HemiTree(indexLocation),...
%             tableResults.R2Cluster4_HemiTree(indexLocation)], 'R2', 'R2AllHemiTree', outputfolder{k}, colorsPerType(k, :))
% 
%         pcaMean = [];
%         pcaSTD = [];
%         for i = 1:size(indexLocationMain,1)
%             spR = split(tableResults.PcaAccuracyMeanSTD_cluster1{indexLocationMain(i)}, ',');
%             pcaMean(i,1) = str2double(spR{1});        
%             pcaSTD(i,1) = str2double(spR{2});
%             spR = split(tableResults.PcaAccuracyMeanSTD_cluster2{indexLocationMain(i)}, ',');
%             pcaMean(i,2) = str2double(spR{1});
%             pcaSTD(i,2) = str2double(spR{2});
%             spR = split(tableResults.PcaAccuracyMeanSTD_cluster3{indexLocationMain(i)}, ',');
%             pcaMean(i,3) = str2double(spR{1});
%             pcaSTD(i,3) = str2double(spR{2});
%             spR = split(tableResults.PcaAccuracyMeanSTD_cluster4{indexLocationMain(i)}, ',');
%             pcaMean(i,4) = str2double(spR{1});
%             pcaSTD(i,4) = str2double(spR{2});
%         end
% 
%         statPlot([pcaMean(:,1),...
%             pcaMean(:,2),...
%             pcaMean(:,3),...
%             pcaMean(:,4)], 'PCA Accuracy', 'PcaAccuracy', outputfolder{k}, colorsPerType(k, :))
% 
%         statPlot([pcaSTD(:,1),...
%             pcaSTD(:,2),...
%             pcaSTD(:,3),...
%             pcaSTD(:,4)], 'PCA STD Accuracy', 'PcaAccuracySTD', outputfolder{k}, colorsPerType(k, :))
% 
%         
%         for j = 1:4
%             f = figure;hold on;
%             shadedErrorBar([], pcaMean(:,j), pcaSTD(:,j), 'lineProps', '-k');
%             plot(1:size(pcaMean,1), tableResults.PcaChanceLevel(indexLocationMain), '--k', 'LineWidth', 2);
%             title(sprintf('Cluster %d, Pca Mean + STD + Chance ', j));
%             xlabel('#Neuron');
%             ylabel('PCA Accuracy');
%             mysave(f, fullfile(outputfolder{k}, sprintf('PCASummary_cluster%d', j)));
%         end
%         
%         f = figure;hold on;
%         title('PCA accuracy Mean and STD')
%         errorbar(mean(pcaMean), mean(pcaSTD));
%         plot(1:size(pcaMean, 2), ones(1, size(pcaMean, 2))*mean(tableResults.PcaChanceLevel(indexLocationMain)), '--k', 'LineWidth', 2)
%         
%         xticks(1:size(pcaMean, 2));
%         xticklabels({'cluster1', 'cluster2', 'cluster3', 'cluster4'});
%         xlim([0,size(pcaMean, 2)+1]);
%         mysave(f, fullfile(outputfolder{k}, 'PCASummary_All'));
%         
%         pcaMeanAll(k) = {pcaMean};
%         pcaSTDAll(k) = {pcaSTD};
%         pcaChance(k) = {mean(tableResults.PcaChanceLevel(indexLocationMain))};
%         
%         text = '';
%         text = strcat(text, sprintf('Cluster1 pass chance %f \\n', mean(pcaMean(:, 1) > tableResults.PcaChanceLevel(indexLocationMain))));
%         text = strcat(text, sprintf('Cluster2 pass chance %f \\n', mean(pcaMean(:, 2) > tableResults.PcaChanceLevel(indexLocationMain))));
%         text = strcat(text, sprintf('Cluster3 pass chance %f \\n', mean(pcaMean(:, 3) > tableResults.PcaChanceLevel(indexLocationMain))));
%         text = strcat(text, sprintf('Cluster4 pass chance %f \\n', mean(pcaMean(:, 4) > tableResults.PcaChanceLevel(indexLocationMain))));
%         
%         fid=fopen(fullfile(outputfolder{k}, 'statisticPCAPassChanceLevel.txt'),'w');
%         fprintf(fid, text);
%         fclose(fid);
% 
%         
%         text = '';
%         text = strcat(text, sprintf('Cluster1 mean %f \\n', mean(tableResults.AnovaTestBetweenAndWithIn_cluster1(indexLocationMain))));
%         text = strcat(text, sprintf('Cluster2 mean %f \\n', mean(tableResults.AnovaTestBetweenAndWithIn_cluster2(indexLocationMain))));
%         text = strcat(text, sprintf('Cluster3 mean %f \\n', mean(tableResults.AnovaTestBetweenAndWithIn_cluster3(indexLocationMain))));
%         text = strcat(text, sprintf('Cluster4 mean %f \\n', mean(tableResults.AnovaTestBetweenAndWithIn_cluster4(indexLocationMain))));
% 
%         fid=fopen(fullfile(outputfolder{k}, 'statisticMeanAnovaSignificantBetweenGroups.txt'),'w');
%         fprintf(fid, text);
%         fclose(fid);
    end
    
%     f = figure;hold on;
%     title('PCA accuracy Mean')
%     groupsBox1 = {};
%     groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster1'};
%     groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster2'};
%     groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster3'};
%     groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-Cluster4'};
%     groupsBox2 = {};
%     groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster1'};
%     groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster2'};
%     groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster3'};
%     groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-Cluster4'};
% 
%     meanAll(1:4) = mean(pcaMeanAll{1});
%     meanAll(5:8) = mean(pcaMeanAll{2});
%     
%     chaAll(1:4) = ones(1,4)*(pcaChance{1});
%     chaAll(5:8) = ones(1,4)*(pcaChance{2});
%     
% 
% 
%     bC = boxchart(categorical(groupsBox1),[pcaMeanAll{1}(:,1); pcaMeanAll{1}(:,2); pcaMeanAll{1}(:,3); pcaMeanAll{1}(:,4)]);
%     
%     bC2 = boxchart(categorical(groupsBox2),[pcaMeanAll{2}(:,1); pcaMeanAll{2}(:,2); pcaMeanAll{2}(:,3); pcaMeanAll{2}(:,4)]);
%     
%   
%     bC.BoxFaceColor = [0,0,0];
%     bC.BoxFaceAlpha = 0.4;
%     bC.MarkerColor = [0,0,0];
%     bC2.BoxFaceColor = [0,0,255] ./ 255;
%     bC2.MarkerColor = [0,0,255] ./ 255;
%     bC2.BoxFaceAlpha = 0.4;
%     
%     plot(meanAll, '-*k');
%     plot(chaAll, '--k')
%     
%     mysave(f, fullfile(outputfolder{1}, 'PCASummary_All_BigAndSmall'));
%     
%     [h1, p1] = ttest2(pcaMeanAll{1}(:,1),  pcaMeanAll{2}(:,1));
%     [h2, p2] = ttest2(pcaMeanAll{1}(:,2),  pcaMeanAll{2}(:,2));
%     [h3, p3] = ttest2(pcaMeanAll{1}(:,3),  pcaMeanAll{2}(:,3));
%     [h4, p4] = ttest2(pcaMeanAll{1}(:,4), pcaMeanAll{2}(:,4));
% 
%     text = '';
%     text = strcat(text, sprintf('Cluster1 h1 %f, p %f \\n', h1,p1));
%     text = strcat(text, sprintf('Cluster2 h1 %f, p %f \\n', h2,p2));
%     text = strcat(text, sprintf('Cluster3 h1 %f, p %f \\n', h3,p3));
%     text = strcat(text, sprintf('Cluster4 h1 %f, p %f \\n', h4,p4));
% 
%     fid=fopen(fullfile(outputfolder{1}, 'statisticTtestPCASmallVsBig.txt'),'w');
%     fprintf(fid, text);
%     fclose(fid);
%     
%     
%     f = figure; hold on;
%     groupsBox1 = {};
%     groupsBox1(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation1))) = {'Big-Cluster1'};
%     groupsBox1(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation1))) = {'Big-Cluster2'};
%     groupsBox1(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation1))) = {'Big-Cluster3'};
%     groupsBox1(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation1))) = {'Big-Cluster4'};
%     groupsBox2 = {};
%     groupsBox2(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation2))) = {'Small-Cluster1'};
%     groupsBox2(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation2))) = {'Small-Cluster2'};
%     groupsBox2(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation2))) = {'Small-Cluster3'};
%     groupsBox2(end+1:end+length(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation2))) = {'Small-Cluster4'};
%     
%     meanAll(1) = mean(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation1));
%     meanAll(2) = mean(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation1));
%     meanAll(3) = mean(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation1));   
%     meanAll(4) = mean(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation1));  
%     meanAll(5) = mean(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation2));
%     meanAll(6) = mean(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation2));    
%     meanAll(7) = mean(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation2));  
%     meanAll(8) = mean(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation2));
%     
% 
%     bC = boxchart(categorical(groupsBox1),[tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation1); ...
%         tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation1);...
%         tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation1);...
%         tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation1)]);
%     
%     bC2 = boxchart(categorical(groupsBox2),[tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation2);...
%         tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation2);...
%         tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation2);...
%         tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation2)]);
%     
%   
%     bC.BoxFaceColor = [0,0,0];
%     bC.BoxFaceAlpha = 0.4;
%     bC.MarkerColor = [0,0,0];
%     bC2.BoxFaceColor = [0,0,255] ./ 255;
%     bC2.MarkerColor = [0,0,255] ./ 255;
%     bC2.BoxFaceAlpha = 0.4;
%     
%     plot(meanAll, '-*k');
%     
%     mysave(f, fullfile(outputfolder{1}, 'R2MaxSmallVsBig.txt'))
%     [h1, p1] = ttest2(tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation1), tableResults.R2MaxBetweenHemiTreeClusre1(indexLocation2));
%     [h2, p2] = ttest2(tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation1), tableResults.R2MaxBetweenHemiTreeClusre2(indexLocation2));
%     [h3, p3] = ttest2(tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation1), tableResults.R2MaxBetweenHemiTreeClusre3(indexLocation2));
%     [h4, p4] = ttest2(tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation1), tableResults.R2MaxBetweenHemiTreeClusre4(indexLocation2));
% 
%     text = '';
%     text = strcat(text, sprintf('Cluster1 h1 %f, p %f \\n', h1,p1));
%     text = strcat(text, sprintf('Cluster2 h1 %f, p %f \\n', h2,p2));
%     text = strcat(text, sprintf('Cluster3 h1 %f, p %f \\n', h3,p3));
%     text = strcat(text, sprintf('Cluster4 h1 %f, p %f \\n', h4,p4));
% 
%     fid=fopen(fullfile(outputfolder{1}, 'statisticTtestR2MaxSmallVsBig.txt'),'w');
%     fprintf(fid, text);
%     fclose(fid);
%     
    calcAnovaForPearsonCompare(matAll, outputfolder{1});
end

function calcAnovaForPearsonCompare(matAll, outputfolder)
    for i = 1:4
        y1 = matAll(matAll(:, 2) == i & matAll(:, 3) == 1& matAll(:, 4) == 1, 1);
        y2 = matAll(matAll(:, 2) == i & matAll(:, 3) == 1& matAll(:, 4) == 2, 1);
        y3 = matAll(matAll(:, 2) == i & matAll(:, 3) == 2& matAll(:, 4) == 1, 1);
        y4 = matAll(matAll(:, 2) == i & matAll(:, 3) == 2& matAll(:, 4) == 2, 1);
        y5 = matAll(matAll(:, 2) == i & matAll(:, 3) == 3& matAll(:, 4) == 1, 1);
        y6 = matAll(matAll(:, 2) == i & matAll(:, 3) == 3& matAll(:, 4) == 2, 1);
        
        sumValues = [];
        sumBW = [];
        sumType = [];
        delta = 60;
        for g = 1:delta:size(y1,1)
            sumValues(end+1) = mean(y1(g:g+10));
            sumBW(end+1) = 1;
            sumType(end+1) = 1;
        end
        
        for g = 1:delta:size(y2,1)
            sumValues(end+1) = mean(y2(g:g+10));
            sumBW(end+1) = 2;
            sumType(end+1) = 1;
        end
        
        for g = 1:delta:size(y3,1)
            sumValues(end+1) = mean(y3(g:g+10));
            sumBW(end+1) = 1;
            sumType(end+1) = 2;
        end     
        
        for g = 1:delta:size(y4,1)
            sumValues(end+1) = mean(y4(g:g+10));
            sumBW(end+1) = 2;
            sumType(end+1) = 2;
        end
        
        for g = 1:delta:size(y5,1)
            sumValues(end+1) = mean(y5(g:g+10));
            sumBW(end+1) = 1;
            sumType(end+1) = 3;
        end     
        
        for g = 1:delta:size(y6,1)
            sumValues(end+1) = mean(y6(g:g+10));
            sumBW(end+1) = 2;
            sumType(end+1) = 3;
        end
        
        [p,tbl,stats] = anovan(sumValues,{sumBW sumType},'model','interaction','varnames',{'B/W','Type'});
        writetable(cell2table(tbl), fullfile(outputfolder, sprintf('AnovaTestByCluster_%d_TOPearsonValuesB_W.csv', i)));
        f = figure; hold on;
        
        meanAll1 = [mean(y1), mean(y2)];
        meanAll2 = [mean(y3), mean(y4)];
        meanAll3 = [mean(y5), mean(y6)];
        stdAll1 = [std(y1), std(y2)];
        stdAll2 = [std(y3), std(y4)];
        stdAll3 = [std(y5), std(y6)];
        
        c1 = cell(1, length(y1));
        c1(:) = {'Big-Within'};
        c2 = cell(1, length(y2));
        c2(:) = {'Big-Between'};
        c3 = cell(1, length(y3));
        c3(:) = {'Small-HZ-Within'};
        c4 = cell(1, length(y4));
        c4(:) = {'Small-HZ-Between'};
        c5 = cell(1, length(y5));
        c5(:) = {'Small-HZTuft-Within'};
        c6 = cell(1, length(y6));
        c6(:) = {'Small-HZTuft-Between'};
        
        bc1 = boxchart(categorical(c1),y1);
        bc2 = boxchart(categorical(c2),y2);
        bc3 = boxchart(categorical(c3),y3);
        bc4 = boxchart(categorical(c4),y4);
        bc5 = boxchart(categorical(c5),y5);
        bc6 = boxchart(categorical(c6),y6);
        
        plot(categorical({'Big-Within', 'Big-Between'}), meanAll1, '-*k');
        plot(categorical({'Small-HZ-Within', 'Small-HZ-Between'}), meanAll2, '-*k');
        plot(categorical({'Small-HZTuft-Within', 'Small-HZTuft-Between'}), meanAll3, '-*k');
        
        mysave(f, fullfile(outputfolder, sprintf('AnovaTestByCluster_%d_TOPearsonValuesB_WFig', i)));
        
        textAnova = '';
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Big = %f, std within Big = %f \\n ',i, mean(y1), std(y1)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Big = %f,  std Between Big = %f \\n ',i, mean(y2), std(y2)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Small HZ = %f, std within Small HZ = %f \\n ',i, mean(y3), std(y3)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Small HZ = %f, std Between Small HZ = %f \\n ',i, mean(y4), std(y4)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean within Small HZTuft = %f, std within Small HZTuft = %f \\n ',i, mean(y5), std(y5)));
        textAnova = strcat(textAnova, sprintf('Cluster %d,  mean Between Small HZTuft = %f, std Between Small HZTuft = %f \\n ',i, mean(y6), std(y6)));
        
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

        [h1,p1] = ttest2(sumValues(sumType==1&sumBW==2), sumValues(sumType==2&sumBW==2), 'Alpha', 0.05);
        [h2,p2] = ttest2(sumValues(sumType==1&sumBW==1),sumValues(sumType==2&sumBW==1), 'Alpha', 0.05);
        [h3,p3] = ttest2(sumValues(sumType==1&sumBW==1),sumValues(sumType==3&sumBW==1), 'Alpha', 0.05);
        [h4,p4] = ttest2(sumValues(sumType==1&sumBW==2),sumValues(sumType==3&sumBW==2), 'Alpha', 0.05);

        tableTtest = table([h1;h2;h3;h4], [p1;p2;p3;p4], 'VariableNames', {'H', 'Pvalue'}, 'RowNames', {'Big_Small_HZ_B','Big_Small_HZ_W','Big_Small_HZTuft_B','Big_Small_HZTuft_W'});
        writetable(tableTtest, fullfile(outputfolder, sprintf('TTestResultsWithinBetweenBigSmallByData%d.csv', i)));

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
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\*\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\*\AVSD*_Depth2_*p*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster1\ByH\AVSD*_Depth2_*p*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster1betweenAll(end+1:end+length(tableBetween.Var2)) = (tableBetween.Var2);
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster1withinAll(end+1:end+length(tablewithin.Var2)) = (tablewithin.Var2);
        end
        
        cluster1withinAll(isnan(cluster1withinAll)) = [];
        cluster1betweenAll(isnan(cluster1betweenAll)) = [];
         
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\*\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\*\AVSD*_Depth2_*p*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster2\ByH\AVSD*_Depth2_*p*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster2betweenAll(end+1:end+length(tableBetween.Var2)) = (tableBetween.Var2);
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster2withinAll(end+1:end+length(tablewithin.Var2)) = (tablewithin.Var2);
        end
        
        cluster2withinAll(isnan(cluster2withinAll)) = [];
        cluster2betweenAll(isnan(cluster2betweenAll)) = [];
         
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\*\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\*\AVSD*_Depth2_*p*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster3\ByH\AVSD*_Depth2_*p*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster3betweenAll(end+1:end+length(tableBetween.Var2)) = (tableBetween.Var2);
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster3withinAll(end+1:end+length(tablewithin.Var2)) = (tablewithin.Var2);
        end
        
        cluster3withinAll(isnan(cluster3withinAll)) = [];
        cluster3betweenAll(isnan(cluster3betweenAll)) = [];
           
        
        cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\*\AVSD*_Depth2_*BetweenMainDepth.csv'];
        cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\*\AVSD*_Depth2_*p*.csv'];
        
        fileBetween = dir(cluster1BetweenFile);
        
        if isempty(fileBetween)
            cluster1BetweenFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\AVSD*_Depth2_*BetweenMainDepth.csv'];
            cluster1WithInFile = [tableResults.RunLocation{indexLocation(i)}, '\cluster4\ByH\AVSD*_Depth2_*p*.csv'];
            fileBetween = dir(cluster1BetweenFile);
        end
        
        tableBetween = readtable(fullfile(fileBetween(1).folder,fileBetween(1).name));
        cluster4betweenAll(end+1:end+length(tableBetween.Var2)) = (tableBetween.Var2);
        
        filewithin = dir(cluster1WithInFile);
        for j = 1:size(filewithin, 1)
            if contains(filewithin(j).folder, '\BetweenSubTrees') || contains(filewithin(j).folder, '\Type2')
                continue;
            end
            tablewithin = readtable(fullfile(filewithin(j).folder,filewithin(j).name));
            cluster4withinAll(end+1:end+length(tablewithin.Var2)) = (tablewithin.Var2);
        end
        
        cluster4withinAll(isnan(cluster4withinAll)) = [];
        cluster4betweenAll(isnan(cluster4betweenAll)) = [];
        
        
    end 
  
    cluster1withinAll_mean = [];
    cluster2withinAll_mean = [];
    cluster3withinAll_mean = [];
    cluster4withinAll_mean = [];
    for g = 1:50:(length(cluster4withinAll))
        cluster4withinAll_mean(end+1) = mean(cluster4withinAll(g:g+10));
        cluster3withinAll_mean(end+1) = mean(cluster3withinAll(g:g+10));
        cluster2withinAll_mean(end+1) = mean(cluster2withinAll(g:g+10));
        cluster1withinAll_mean(end+1) = mean(cluster1withinAll(g:g+10));
    end
    
    
    cluster1BetweenAll_mean = [];
    cluster2BetweenAll_mean = [];
    cluster3BetweenAll_mean = [];
    cluster4BetweenAll_mean = [];
    for g = 1:50:(length(cluster4betweenAll))
        cluster4BetweenAll_mean(end+1) = mean(cluster4betweenAll(g:g+10));
        cluster3BetweenAll_mean(end+1) = mean(cluster3betweenAll(g:g+10));
        cluster2BetweenAll_mean(end+1) = mean(cluster2betweenAll(g:g+10));
        cluster1BetweenAll_mean(end+1) = mean(cluster1betweenAll(g:g+10));
    end
    
    [h1,p1] = ttest2(cluster1withinAll_mean, cluster1BetweenAll_mean, 'Alpha', alphValue);
    [h2,p2] = ttest2(cluster2withinAll_mean, cluster2BetweenAll_mean, 'Alpha', alphValue);
    [h3,p3] = ttest2(cluster3withinAll_mean, cluster3BetweenAll_mean, 'Alpha', alphValue);   
    [h4,p4] = ttest2(cluster4BetweenAll_mean,cluster4withinAll_mean, 'Alpha', alphValue);
       
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