function RunSummaryForHemiTreePearson_AllAndExcludedEvents_191021()
    mantelRunnerLocation = '\\jackie-analysis\e\Shay\RunnersLocationSummary.xlsx';
    sheetName = 'RunOnlyTuft';
    outputfolder(1) = {'\\jackie-analysis\e\Shay\StatisticSummary\HemiTreePostSummary\Big_19_10_Final1\'};
    mkdir(outputfolder{1});
    classification(1) = 1;
    
    outputfolder(2) = {'\\jackie-analysis\e\Shay\StatisticSummary\HemiTreePostSummary\Small_19_10_Final1\'};
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
    
    r2All = {};
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
       
        [~, bslop, r2] = BetweenVsWithin(indexLocationMain,indexLocation, tableResults, outputfolder{k}, 0.05);
        
        slopAll(k) = {bslop};
        r2All(k) = {r2};
       
        pcaMean = [];
        pcaSTD = [];
        for i = 1:size(indexLocationMain,1)
            spR = split(tableResults.PcaAccuracyMeanSTD_All{indexLocationMain(i)}, ',');
            pcaMean(i,1) = str2double(spR{1});        
            pcaSTD(i,1) = str2double(spR{2}); 
        end

        statPlot([pcaMean(:,1),pcaMean(:,1)], 'PCA Accuracy', 'PcaAccuracy', outputfolder{k}, colorsPerType(k, :))

        statPlot([pcaSTD(:,1),pcaSTD(:,1)], 'PCA STD Accuracy', 'PcaAccuracySTD', outputfolder{k}, colorsPerType(k, :))

        
        for j = 1
            f = figure;hold on;
            shadedErrorBar([], pcaMean(:,j), pcaSTD(:,j), 'lineProps', '-k');
            plot(1:size(pcaMean,1), tableResults.PcaChanceLevel(indexLocationMain), '--k', 'LineWidth', 2);
            title(sprintf('Cluster All, Pca Mean + STD + Chance '));
            xlabel('#Neuron');
            ylabel('PCA Accuracy');
            mysave(f, fullfile(outputfolder{k}, sprintf('PCASummary_clusterAll')));
        end
        
        f = figure;hold on;
        title('PCA accuracy Mean and STD')
        errorbar(mean(pcaMean), mean(pcaSTD));
        plot(mean(pcaMean), '--*');
        plot(0:2, ones(1, size(pcaMean, 2))*mean(tableResults.PcaChanceLevel(indexLocationMain)), '--k', 'LineWidth', 2)
        
        xticks([1]);
        xticklabels({'All'});
        xlim([0,2]);
        mysave(f, fullfile(outputfolder{k}, 'PCASummary_All'));
        
        pcaMeanAll(k) = {pcaMean};
        pcaSTDAll(k) = {pcaSTD};
        pcaChance(k) = {mean(tableResults.PcaChanceLevel(indexLocationMain))};
        
        text = '';
        text = strcat(text, sprintf('All pass chance %f \\n', mean(pcaMean(:, 1) > tableResults.PcaChanceLevel(indexLocationMain))));
        
        fid=fopen(fullfile(outputfolder{k}, 'statisticPCAPassChanceLevel.txt'),'w');
        fprintf(fid, text);
        fclose(fid);
    end
    
    f = figure;hold on;
    title('Slop')
    groupsBox1 = {};
    groupsBox1(end+1:end+size(slopAll{1}, 1)) = {'Big-All'};
    groupsBox1(end+1:end+size(slopAll{1}, 1)) = {'Big-EB'};
    groupsBox2 = {};
    groupsBox2(end+1:end+size(slopAll{2}, 1)) = {'Small-All'};
    groupsBox2(end+1:end+size(slopAll{2}, 1)) = {'Small-EB'};

    bC = boxchart(categorical(groupsBox1),[slopAll{1}(:,1); slopAll{1}(:,2)]);    
    bC2 = boxchart(categorical(groupsBox2),[slopAll{2}(:,1); slopAll{2}(:,2)]);
    
    bC.BoxFaceColor = [0,0,0];
    bC.BoxFaceAlpha = 0.4;
    bC.MarkerColor = [0,0,0];
    bC2.BoxFaceColor = [0,0,255] ./ 255;
    bC2.MarkerColor = [0,0,255] ./ 255;
    bC2.BoxFaceAlpha = 0.4;
    
    plot([mean(slopAll{1}),mean(slopAll{2})], '-*k');
        
    mysave(f, fullfile(outputfolder{1}, 'SlopSummary_All_BigAndSmall'));
    
    [p1, h1] = ranksum(slopAll{1}(:,1),  slopAll{2}(:,1));
    [p2, h2] = ranksum(slopAll{1}(:,2),  slopAll{2}(:,2));

    text = '';
    text = strcat(text, sprintf('All h1 %f, p %f \\n', h1,p1));
    text = strcat(text, sprintf('EB h1 %f, p %f \\n', h2,p2));

    fid=fopen(fullfile(outputfolder{1}, 'statisticranksumSlopSmallVsBig.txt'),'w');
    fprintf(fid, text);
    fclose(fid);
    
    f = figure;hold on;
    title('R2')
    groupsBox1 = {};
    groupsBox1(end+1:end+size(r2All{1}, 1)) = {'Big-All'};
    groupsBox1(end+1:end+size(r2All{1}, 1)) = {'Big-EB'};
    groupsBox2 = {};
    groupsBox2(end+1:end+size(r2All{2}, 1)) = {'Small-All'};
    groupsBox2(end+1:end+size(r2All{2}, 1)) = {'Small-EB'};

    bC = boxchart(categorical(groupsBox1),[r2All{1}(:,1); r2All{1}(:,2)]);    
    bC2 = boxchart(categorical(groupsBox2),[r2All{2}(:,1); r2All{2}(:,2)]);
    
    bC.BoxFaceColor = [0,0,0];
    bC.BoxFaceAlpha = 0.4;
    bC.MarkerColor = [0,0,0];
    bC2.BoxFaceColor = [0,0,255] ./ 255;
    bC2.MarkerColor = [0,0,255] ./ 255;
    bC2.BoxFaceAlpha = 0.4;
    
    plot([mean(r2All{1}),mean(r2All{2})], '-*k');
    
    mysave(f, fullfile(outputfolder{1}, 'R2Summary_All_BigAndSmall'));
    
    [p1, h1] = ranksum(r2All{1}(:,1),  r2All{2}(:,1));
    [p2, h2] = ranksum(r2All{1}(:,2),  r2All{2}(:,2));

    text = '';
    text = strcat(text, sprintf('All h1 %f, p %f \\n', h1,p1));
    text = strcat(text, sprintf('EB h1 %f, p %f \\n', h2,p2));

    fid=fopen(fullfile(outputfolder{1}, 'statisticranksumR2SmallVsBig.txt'),'w');
    fprintf(fid, text);
    fclose(fid);
    
    f = figure;hold on;
    title('PCA accuracy Mean')
    groupsBox1 = {};
    groupsBox1(end+1:end+size(pcaMeanAll{1}, 1)) = {'Big-All'};
    groupsBox2 = {};
    groupsBox2(end+1:end+size(pcaMeanAll{2}, 1)) = {'Small-All'};

    meanAll(1) = mean(pcaMeanAll{1});
    meanAll(2) = mean(pcaMeanAll{2});
    
    chaAll(1) = ones(1,1)*(pcaChance{1});
    chaAll(2) = ones(1,1)*(pcaChance{2});
    


    bC = boxchart(categorical(groupsBox1),[pcaMeanAll{1}(:,1)]);
    
    bC2 = boxchart(categorical(groupsBox2),[pcaMeanAll{2}(:,1)]);
    
  
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

    text = '';
    text = strcat(text, sprintf('All h1 %f, p %f \\n', h1,p1));

    fid=fopen(fullfile(outputfolder{1}, 'statisticranksumPCASmallVsBig.txt'),'w');
    fprintf(fid, text);
    fclose(fid);
end

function [matAll, bslop, r2] = BetweenVsWithin(indexLocation,indexHemi, tableResults, outputfolder, alphValue)
    clusterAllbetweenAll = [];
    clusterAllwithinAll = [];
    clusterExcludebigbetweenAll = [];
    clusterExcludebigwithinAll = [];

    tableResults.b1_clusterAll = nan(size(tableResults, 1),1);
    tableResults.b1_clusterExcludebig = nan(size(tableResults, 1),1);   
    tableResults.R2_clusterAll = nan(size(tableResults, 1),1);
    tableResults.R2_clusterExcludebig = nan(size(tableResults, 1),1);
      
    for i = 1:size(indexLocation, 1)
        clusterAllFile = [tableResults.RunLocation{indexLocation(i)}, '\All\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*.csv'];
        
        fileList = dir(clusterAllFile);
        
        xVals = [];
        yVals = [];
        x1Vals = {};
        y1Vals = {};
        in_k = 1;
        
        temp1 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
            
            if contains(fileList(j).name, 'BetweenMainDepth')
                clusterAllbetweenAll(end+1) = nanmean(tableR.Var2);
            else
                temp1(end+1:end+length(tableR.Var2)) = tableR.Var2;
                x1Vals(in_k) = {tableR.Var1};
                y1Vals(in_k) = {tableR.Var2};
                 in_k = in_k + 1;
            end
            
            
            xVals(end+1:end+length(tableR.Var1)) = tableR.Var1;
            yVals(end+1:end+length(tableR.Var2)) = tableR.Var2;
        end
        clusterAllwithinAll(end+1) = nanmean(temp1);
            
        clusterAllwithinAll(isnan(clusterAllwithinAll)) = [];
        clusterAllbetweenAll(isnan(clusterAllbetweenAll)) = [];
         
        mdAll1 = fitglm(xVals ./ max(xVals), yVals);
        b1_clusterAll = mdAll1.Coefficients.Estimate(2);
        R2_clusterAll = mdAll1.Rsquared.Ordinary;
         
        tableResults.b1_clusterAll(indexLocation(i)) = b1_clusterAll;
        tableResults.R2_clusterAll(indexLocation(i)) = R2_clusterAll;
        
        clusterExcludeBigFile = [tableResults.RunLocation{indexLocation(i)}, '\All_ExcludeBigEvents\ByH\BetweenAndWithinSubTrees\AVSD*_Depth2_*.csv'];
        
        fileList = dir(clusterExcludeBigFile);
        
        xVals = [];
        yVals = [];
       
        x1Vals = {};
        y1Vals = {};
        in_k = 1;
        
        temp2 = [];
        for j = 1:size(fileList, 1)
            tableR = readtable(fullfile(fileList(j).folder,fileList(j).name));
            
            if contains(fileList(j).name, 'BetweenMainDepth')
                clusterExcludebigbetweenAll(end+1) = nanmean(tableR.Var2);
            else
                temp2(end+1:end+length(tableR.Var2)) = tableR.Var2;
                x1Vals(in_k) = {tableR.Var1};
                y1Vals(in_k) = {tableR.Var2};
                in_k = in_k + 1;
            
            end
            
            
            xVals(end+1:end+length(tableR.Var1)) = tableR.Var1;
            yVals(end+1:end+length(tableR.Var2)) = tableR.Var2;
        end
        clusterExcludebigwithinAll(end+1) = nanmean(temp2);
            
        clusterExcludebigwithinAll(isnan(clusterExcludebigwithinAll)) = [];
        clusterExcludebigbetweenAll(isnan(clusterExcludebigbetweenAll)) = [];
         
        
        mdAll1 = fitglm(xVals ./ max(xVals), yVals);
        b1_clusterExcludebig = mdAll1.Coefficients.Estimate(2);
        R2_clusterExcludebig = mdAll1.Rsquared.Ordinary;
        
        tableResults.b1_clusterExcludebig(indexLocation(i)) = b1_clusterExcludebig;
        tableResults.R2_clusterExcludebig(indexLocation(i)) = R2_clusterExcludebig;
    end
    
    save(fullfile(outputfolder, 'tablewithB1'), 'tableResults')   
    
    matResultsV = [clusterAllwithinAll';clusterAllbetweenAll';clusterExcludebigwithinAll';clusterExcludebigbetweenAll'];
    
    matResultsCluster = [ones(length(clusterAllwithinAll),1); ones(length(clusterAllbetweenAll),1);...
        ones(length(clusterExcludebigwithinAll),1)*2;ones(length(clusterExcludebigbetweenAll),1)*2];
    
    matResultsType = [ones(length(clusterAllwithinAll),1); ones(length(clusterAllbetweenAll),1);...
        ones(length(clusterExcludebigwithinAll),1);ones(length(clusterExcludebigbetweenAll),1)];
    
    matResultsBW = [ones(length(clusterAllwithinAll),1); ones(length(clusterAllbetweenAll),1)*2;...
        ones(length(clusterExcludebigwithinAll),1);ones(length(clusterExcludebigbetweenAll),1)*2];
    
    matAll = [matResultsV, matResultsCluster, matResultsType, matResultsBW];
    
    f = figure; hold on;
    title('Linear regression slop (b1)');
    groupsBox1 = {};
    groupsBox1(end+1:end+length(indexLocation)) = {'All'};
    groupsBox1(end+1:end+length(indexLocation)) = {'Excludebig'};
    meanAll = [];
    bC = boxchart(categorical(groupsBox1), [tableResults.b1_clusterAll(indexLocation); ...
        tableResults.b1_clusterExcludebig(indexLocation)]);
    meanAll(1) = mean(tableResults.b1_clusterAll(indexLocation));
    meanAll(2) = mean(tableResults.b1_clusterExcludebig(indexLocation));  
    plot(meanAll, '-*k');
    ylim([-1.2, 0.2]);
    mysave(f, fullfile(outputfolder, 'B1_slop'));   
    
    statPlot([tableResults.b1_clusterAll(indexLocation), ...
        tableResults.b1_clusterExcludebig(indexLocation)], 'slop(b1)', 'Linear regression slop (b1)', outputfolder, [0,0,128]./255);
    
    statPlot([tableResults.R2_clusterAll(indexLocation), ...
        tableResults.R2_clusterExcludebig(indexLocation)], 'R2', 'Linear regression R2', outputfolder, [0,0,128]./255);
    
    
    bslop = [tableResults.b1_clusterAll(indexLocation), ...
        tableResults.b1_clusterExcludebig(indexLocation)];   
    r2 = [tableResults.R2_clusterAll(indexLocation), ...
        tableResults.R2_clusterExcludebig(indexLocation)];
end

function statPlot(arrayT, ylabelName, TitleV, outputfolder, colorType)
    f = figure; hold on;
    b = boxchart(arrayT);
    ylim([0,1]);
    xticklabels({'All', 'EB'});
    b.BoxFaceColor = colorType;
    b.BoxFaceAlpha = 0.4;
    b.MarkerColor = colorType;
    ylabel(ylabelName);
    plot(mean(arrayT), '-*k');
    mysave(f, fullfile(outputfolder, TitleV));
    mResult = mean(arrayT);
    stdR = std(arrayT);
    textAnova = '';

    textAnova = strcat(textAnova, sprintf('All mean: %.4f, std: %.4f \\n ', mResult(1), stdR(1)));
    textAnova = strcat(textAnova, sprintf('EB mean: %.4f, std: %.4f \\n ', mResult(2), stdR(2)));            

    fid=fopen(fullfile(outputfolder, ['summaryMandstd_', TitleV,'.txt']),'w');
    fprintf(fid, textAnova);
    fclose(fid);
    
    textAnova = '';
    [p,~,statsM] = anova1(arrayT, {'All', 'EB'}, 'off');
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