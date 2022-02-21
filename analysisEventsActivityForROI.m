function analysisEventsActivityForROI(gROI, allEventsTable, selectedROI, selectedROISplitDepth, outputpath, fileName, clusterCount, isMainDepth, mainDepthCount)
    outputpath = [outputpath, '\PrecentagePerSubTrees\'];
    mkdir(outputpath);
    
    classes = unique(selectedROISplitDepth);
    classes(classes == -1) = [];
    fileID = fopen([outputpath, '\' fileName 'Sum_behaveAnalysis.txt'],'w');
    fileIDSec = fopen([outputpath, '\' fileName 'Precentage_behaveAnalysis.txt'],'w');
    
    for i_cluster = 0:clusterCount
        sumActivityForROI_noW = zeros(1, length(selectedROI));
        sumActivityForROI_W_PksH = zeros(1, length(selectedROI));
        sumActivityForROI_W_Precentage = zeros(1, length(selectedROI));

        countEvents_noW = 0;
        countEvents_P = 0;
        
        for i_e = 1:size(allEventsTable, 1)
            if allEventsTable.clusterByH(i_e) == i_cluster || i_cluster == 0
                sumActivityForROI_W_PksH(allEventsTable.roisEvent{i_e} == 1) = sumActivityForROI_W_PksH(allEventsTable.roisEvent{i_e} == 1) + 1*allEventsTable.H(i_e);
                sumActivityForROI_noW(allEventsTable.roisEvent{i_e} == 1) = sumActivityForROI_noW(allEventsTable.roisEvent{i_e} == 1) + 1;        
                countEvents_noW = countEvents_noW + 1;
            end
            
            if allEventsTable.clusterByRoiPrecantage(i_e) == i_cluster || i_cluster == 0
                sumActivityForROI_W_Precentage(allEventsTable.roisEvent{i_e} == 1) = sumActivityForROI_W_Precentage(allEventsTable.roisEvent{i_e} == 1) + 1*allEventsTable.roiPrecantage(i_e);
                countEvents_P = countEvents_P + 1;
            end
        end
        
        plotResultsA(outputpath, [fileName, '_noWeight', '_sum'], i_cluster, sumActivityForROI_noW, classes, selectedROISplitDepth, 'no weight, Average Sum', gROI, fileID,  isMainDepth, mainDepthCount);
        plotResultsA(outputpath, [fileName, '_WeightHigh', '_sum'], i_cluster, sumActivityForROI_W_PksH, classes, selectedROISplitDepth, 'weight by High(Pks), Average Sum', gROI, fileID,  isMainDepth, mainDepthCount);
        plotResultsA(outputpath, [fileName, '_WeightPrecentage', '_sum'], i_cluster, sumActivityForROI_W_Precentage, classes, selectedROISplitDepth, 'weight by Precentage, Average Sum', gROI, fileID,  isMainDepth, mainDepthCount);   
        
        if countEvents_noW ~= 0 
            sumActivityForROI_noW = sumActivityForROI_noW ./ countEvents_noW;
            sumActivityForROI_W_PksH = sumActivityForROI_W_PksH ./ countEvents_noW;
            
            plotResultsA(outputpath, [fileName, '_noWeight', '_precentage'], i_cluster, sumActivityForROI_noW, classes, selectedROISplitDepth, 'no weight, Average Precentage', gROI, fileIDSec,  isMainDepth, mainDepthCount);
            plotResultsA(outputpath, [fileName, '_WeightHigh', '_precentage'], i_cluster, sumActivityForROI_W_PksH, classes, selectedROISplitDepth, 'weight by High(Pks), Average Precentage', gROI, fileIDSec,  isMainDepth, mainDepthCount);
        
        end
        
        if countEvents_P ~= 0
            sumActivityForROI_W_Precentage = sumActivityForROI_W_Precentage ./ countEvents_P;
            plotResultsA(outputpath, [fileName, '_WeightPrecentage', '_precentage'], i_cluster, sumActivityForROI_W_Precentage, classes, selectedROISplitDepth, 'weight by Precentage, Average Precentage', gROI, fileIDSec,  isMainDepth, mainDepthCount);   
        end
        
    end
    
    fclose(fileID);
end

function plotResultsA(outputpath, fileName, i_cluster, sumActivityForROI, classes, selectedROISplitDepth, titleAdd, gROI, fileID,  isMainDepth, mainDepthCount)
    formatSpec = 'cluster %d, type: %s, name: %s, mean: %d, std: %d, roi#: %d\n';
    
    fig = figure;
    hold on;

    %         all ROI        
    meanAll = mean(sumActivityForROI);
    stdAll = std(sumActivityForROI);
    x_tick{1} = 'allRoi';

    errorbar(1, meanAll,stdAll, 'o', 'Color', 'black', 'MarkerSize', 6, 'MarkerEdgeColor','black','MarkerFaceColor','black');
    text(1,meanAll,sprintf(' mean: %.2f,\n std: %.2f\n roi#: %d', meanAll, stdAll, length(sumActivityForROI)), 'FontSize', 8);
    indexT = 3;
    
    fprintf(fileID, formatSpec,i_cluster,titleAdd,'allROI', meanAll, stdAll, length(sumActivityForROI));
    anovaY = [];
    anovaGr = {};
    
    for s_t = 1:length(classes)
        tmp(s_t) = {sumActivityForROI(classes(s_t) == selectedROISplitDepth)};
        sum_c = sum(classes(s_t) == selectedROISplitDepth);
        m_c = mean(sumActivityForROI(classes(s_t) == selectedROISplitDepth));
        std_c = std(sumActivityForROI(classes(s_t) == selectedROISplitDepth));
        
        color = getTreeColor('within', s_t,  isMainDepth, mainDepthCount);
        name(s_t) = gROI.Nodes.Name(classes(s_t));
        x_tick{end + 1} = name{s_t};
        errorbar(indexT, m_c,std_c,'o', 'Color', color, 'MarkerSize', 6, 'MarkerEdgeColor',color,'MarkerFaceColor',color);
        text(indexT, m_c,sprintf(' mean: %.2f,\n std: %.2f\n roi#: %d', m_c, std_c, sum_c), 'FontSize', 8);
        indexT = indexT + 2;
        
        anovaY((end+1):(end+ size(tmp{s_t}, 2))) = tmp{s_t};
        anovaGr((end+1):(end+ size(tmp{s_t}, 2))) = name(s_t);
        
        
        fprintf(fileID, formatSpec,i_cluster,titleAdd,name{s_t}, m_c, std_c, sum_c);
    end

    title({['cluster: ' num2str(i_cluster) ' , Type: ' titleAdd]});
    xticks([1:2:(indexT - 2)]);
    xticklabels(x_tick);
    xtickangle(90);
    xlim([0, (indexT)]);
    ylabel('average number of events');

    mysave(fig, [outputpath, '\' fileName '_cluster_' num2str(i_cluster) 'analysisR']);

%     --------------------------------------------------------------------------------------    
    statisticTestResults = array2table(cell(0,8));
    statisticTestResults.Properties.VariableNames = {'TestType', 'GroupName1', 'GroupName2', 'H_value', 'P_value', 'CI_Upper', 'CI_Low', 'MeanDiff'};
    counterStat = 1;
    
    for i = 1:length(tmp)
        for k = (i+1):length(tmp)
            [h,p] = ttest2(tmp{i}, tmp{k});
            statisticTestResults.TestType(counterStat) = {'Ttest'};
            statisticTestResults.GroupName1(counterStat) = name(i);
            statisticTestResults.GroupName2(counterStat) = name(k);
            statisticTestResults.H_value(counterStat) = {h};
            statisticTestResults.P_value(counterStat) = {p};
            counterStat =counterStat + 1; 
        end
    end
    
    
    [~,~,anova_stats] = anova1(anovaY, anovaGr);
    
    if length(classes) > 1 && anova_stats.df ~= 0
        [c,~,h,gnames] = multcompare(anova_stats);

        for c_ind = 1:size(c, 1)
            statisticTestResults.TestType(counterStat) = {'One way Anova'};
            statisticTestResults.GroupName1(counterStat) = gnames(c(c_ind, 1));
            statisticTestResults.GroupName2(counterStat) = gnames(c(c_ind, 2));
            statisticTestResults.CI_Upper(counterStat) = {c(c_ind, 5)};
            statisticTestResults.CI_Low(counterStat) =  {c(c_ind, 3)};
            statisticTestResults.MeanDiff(counterStat) = {c(c_ind, 4)};
            statisticTestResults.P_value(counterStat) = {c(c_ind, 6)};
            counterStat = counterStat + 1;
        end

        fileName5 = [outputpath, '\AnovaPlotCluster', num2str(i_cluster), '_', fileName];
        mysave(h, fileName5);
    end
    
    writetable(statisticTestResults, [outputpath, '\StatisticsCluster_', num2str(i_cluster), '_', fileName, '.csv'])

end