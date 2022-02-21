function plotAndCalcStatisticTest(classesColorName, classesColor, resultsT, leg, legColor, outputpath, classesM, clusterType, roiActivityPeakSize, distType)
    statisticTestResults = array2table(cell(0,12));
    statisticTestResults.Properties.VariableNames = {'Welch_Fval', 'Welch_Pval', 'Welch_df1',...
        'Welch_df2', 'TestType', 'GroupName1', 'GroupName2', 'H_value', 'P_value',...
        'CI_Upper', 'CI_Low', 'MeanDiff'};
     
    counterStat = 1;
    
    %     Plot ROI Activity VS Tree Distance
    fig3 = figure;
    hold on;

    title({'ROI Activity', distType});
    legend(leg, legColor);

    matR = [];
    groupR = [];
    labelsR = [];
    labelsSummary = {'Group', 'Mean', 'Std', 'maxPval', 'minPval'};
    colorsG = [];
    anovaY = [];
    anovaGr = {};
    indexGroup = 1;
    for clr1 = 1:length(classesColorName)
        tempR = resultsT.(classesColorName{clr1});

        if isempty(tempR)
            continue;
        end

        anovaY((end+1):(end+ size(tempR, 1))) = tempR(:, 2);
        anovaGr((end+1):(end+ size(tempR, 1))) = classesColorName(clr1);
        
        matR = [matR; tempR(:,2)];
        groupR = [groupR; ones(size(tempR(:,2)))*indexGroup];

        indexGroup = indexGroup + 1;

        labelsR{end + 1} = classesColorName{clr1};
      
        temp2 = tempR(:,2);
        temp2(isnan(temp2)) = [];
        labelsSummary(end+1, 1) = classesColorName(clr1);
        labelsSummary(end, 2) = {mean(temp2)};
        labelsSummary(end, 3) = {std(temp2)};
        labelsSummary(end, 4) = {max(tempR(:,3))};
        labelsSummary(end, 5) = {min(tempR(:,3))};
        
        colorsG(end+1, :) = classesColor{clr1};
    end
    
    boxplot(matR, groupR, 'Labels', labelsR, 'Colors', colorsG);        
    xtickangle(90);   
    ylabel(['Roi Activity ', distType]);
    fileName4 = [outputpath, '\BOXPlot_numofTreeDepth', num2str(length(classesM)), roiActivityPeakSize, '_', num2str(clusterType)];
    mysave(fig3, fileName4);
    
    matR2 = matR;
    matR2(isnan(matR2)) = [];
    labelsSummary(end+1, 1) = {'all'};
    labelsSummary(end, 2) = {mean(matR2, 'all')};
    labelsSummary(end, 3) = {std(matR2, 0, 'all')};
    labelsSummary(end, 4) = {nan};
    labelsSummary(end, 5) = {nan};

    statisticRegResults = array2table(labelsSummary);
    writetable(statisticRegResults, [outputpath, '\Summary_numofTreeDepth', num2str(length(classesM)), roiActivityPeakSize, '_', num2str(clusterType), '.csv'])
    
    errorBarSummary = figure;
    colorsG(end+1, :) = [0,0,0];
    hold on;
    for rIndex = 2:size(labelsSummary, 1)
        isSig = '';
        if labelsSummary{rIndex, 4} <= 0.01
            isSig = '**';
        elseif labelsSummary{rIndex, 4} <= 0.05
            isSig = '*';
        end
        
        errorbar(rIndex, labelsSummary{rIndex, 2},labelsSummary{rIndex, 3},'o', 'Color', colorsG(rIndex-1, :), 'MarkerSize', 6, 'MarkerEdgeColor',colorsG(rIndex-1, :),'MarkerFaceColor',colorsG(rIndex-1, :));
        text(rIndex, labelsSummary{rIndex, 2},sprintf(' mean: %.2f,\n std: %.2f\n %s',...
            labelsSummary{rIndex, 2}, labelsSummary{rIndex, 3},isSig), 'FontSize', 8);
    end
    
    
    title({'Activity Correlation Summary'});
    xticklabels(labelsSummary(2:end, 1));
    xtickangle(90);
    xticks([2:(size(labelsSummary, 1))]);
    xlim([0, length(labelsSummary)+1])
    mysave(errorBarSummary, [outputpath, '\MeanSummaryPlot_', num2str(length(classesM))]);
 
    
%     Check Welch Test
    [p, F, df1, df2] = wanova(anovaY', groupR);
    statisticTestResults.Welch_Fval(counterStat) = {F};
    statisticTestResults.Welch_Pval(counterStat) = {p};
    statisticTestResults.Welch_df1(counterStat) = {df1};
    statisticTestResults.Welch_df2(counterStat) = {df2};
    counterStat = counterStat + 1;
    if p < 0.05
        %  Ttest Results
        formatSpec = '\n SubTrees compared: %s vs %s, ttest results: \n h = %d, pValue = %d \n\r\n ';
        ttest_results_str = '';
        ttest_results_str = strcat(ttest_results_str, 'The Welch Test pass significantly, So we can compare groups');
        for clr1 = 1:length(classesColorName)
            temp_1 = resultsT.(classesColorName{clr1})(:,2);
            for clr2 = (clr1+1):length(classesColorName)
                temp_2 = resultsT.(classesColorName{clr2})(:,2);

                [h,p] = ttest2(temp_1, temp_2, 'Vartype', 'unequal');
                ttest_results_str = strcat(ttest_results_str,...
                    sprintf(formatSpec, classesColorName{clr1}, classesColorName{clr2}, h, p));

                statisticTestResults.TestType(counterStat) = {'Ttest'};
                statisticTestResults.GroupName1(counterStat) = classesColorName(clr1);
                statisticTestResults.GroupName2(counterStat) = classesColorName(clr2);
                statisticTestResults.H_value(counterStat) = {h};
                statisticTestResults.P_value(counterStat) = {p};
                counterStat =counterStat + 1;
            end
        end
    else
        ttest_results_str = 'The Welch Test did not pass, So we can not compare groups';
    end
            
    
%     --------------------------------------------------------------------------

    import mlreportgen.ppt.*
    ppt = Presentation([outputpath '\WelchResults_numofTreeDepth', num2str(length(classesM)),  roiActivityPeakSize, '_', num2str(clusterType)],...
        'AnalysisTtest.potm');
    open(ppt);
    currentResultsSlide= add(ppt, 'TtestLyout');
    
    replace(currentResultsSlide.Children(2), Picture([fileName4 '.tif']));       
    replace(currentResultsSlide.Children(1), Paragraph(ttest_results_str)); 
    replace(currentResultsSlide.Children(3), Paragraph(['Ttest Results ', roiActivityPeakSize, '_', num2str(clusterType)]));       
     
% %     ------------------------------------------------------------------------------
% %  one way anova
% 
    if (length(labelsR) > 1)
        [pA,tb1A,anova_stats] = anova1(anovaY, anovaGr);
        if anova_stats.df ~= 0
            ttest_results_str = '';
            
            if pA < 0.05
                ttest_results_str = strcat(ttest_results_str, 'The Anova Test pass significantly, So we can compare groups');
                [c,~,h,gnames] = multcompare(anova_stats);
            
                formatSpecAnova = '\n SubTrees compared: %s vs %s, low_CI = %d mean diff = %d upper_CI = %d pValue = %d \n\r\n ';
                for c_ind = 1:size(c, 1)
                    ttest_results_str = strcat(ttest_results_str,...
                            sprintf(formatSpecAnova, gnames{c(c_ind, 1)}, gnames{c(c_ind, 2)}, c(c_ind, 3), c(c_ind, 4), c(c_ind, 5), c(c_ind, 6)));       

                    statisticTestResults.TestType(counterStat) = {'One way Anova'};
                    statisticTestResults.GroupName1(counterStat) = gnames(c(c_ind, 1));
                    statisticTestResults.GroupName2(counterStat) = gnames(c(c_ind, 2));
                    statisticTestResults.CI_Upper(counterStat) = {c(c_ind, 5)};
                    statisticTestResults.CI_Low(counterStat) =  {c(c_ind, 3)};
                    statisticTestResults.MeanDiff(counterStat) = {c(c_ind, 4)};
                    statisticTestResults.P_value(counterStat) = {c(c_ind, 6)};
                    counterStat = counterStat + 1;
                    
                    fileName5 = [outputpath, '\AnovaPlot_numofTreeDepth', num2str(length(classesM)), roiActivityPeakSize, '_', num2str(clusterType)];
                    mysave(h, fileName5);
                    replace(currentResultsSlide.Children(2), Picture([fileName5 '.tif']));       
            
                end
            else
                ttest_results_str = strcat(ttest_results_str, 'The Anova Test did not pass , So we can not compare groups');
            end
            
            currentResultsSlide= add(ppt, 'AnovaLyout');

            replace(currentResultsSlide.Children(1), Paragraph(ttest_results_str)); 
            replace(currentResultsSlide.Children(3), Paragraph(['Anova one way Results ',roiActivityPeakSize, '_', num2str(clusterType)])); 
            replace(currentResultsSlide.Children(4), Picture([fileName4 '.tif']));       
        end
    end
    
    close(ppt);
   
    writetable(statisticTestResults, [outputpath, '\Statistics_numofTreeDepth', num2str(length(classesM)), roiActivityPeakSize, '_', num2str(clusterType), '.csv'])
end