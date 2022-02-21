function GLMPostAnalysisSummary()
    glmLocationSummary = '\\jackie-analysis10\e\Shay\RunnersLocationSummary.xlsx';
    sheetName = 'RunGLM';
    timeSeg = 1;
    
    classification = 2;
    
    
    passThreshold = 0.15;
    tableResultsAll = readtable(glmLocationSummary,'Sheet',sheetName);
    
    indexIncludeSummary =  tableResultsAll.includeSummaryAll == 1 & tableResultsAll.classification == classification;
    tableResults = tableResultsAll(indexIncludeSummary,:);
        
    outputpathAll = '\\jackie-analysis\e\Shay\StatisticSummary\GLMSummary_Final_1402_FinalNormal';
    
    outputpath = [outputpathAll, '\', num2str(classification)];
    
    for i = 1:size(tableResults, 1)
        resultsByR2.(tableResults.Type{i}).r = [];
        resultsByAnimal.(sprintf('Animal%s', tableResults.Animal{i})).(tableResults.Neuron{i}).(tableResults.Type{i}).sum = [];
        
%         runData = load([tableResults.RunLocation{i}, '\glmResHist.mat']);
%         for z = 1:length(runData.eventsTypes)
%             resultsByCont.(tableResults.Type{i}).(runData.eventsTypes{z}).v = 0;
%             resultsByCont.(tableResults.Type{i}).(runData.eventsTypes{z}).c = 0;
%         end
    end
    
    for i = 1:size(tableResults, 1)
        totalROIs = 0;
        passROIs = 0;
        
        if isfile([tableResults.RunLocation{i}, '\glmResHist.mat'])
            runData = load([tableResults.RunLocation{i}, '\glmResHist.mat']);
        elseif isfile([tableResults.RunLocation{i}, '\glmres.mat'])
            runData = load([tableResults.RunLocation{i}, '\glmres.mat']);
        else
            continue;
        end
        totalROIs = totalROIs + size(runData.R2full_te{timeSeg}, 1);
        
        for k=1:(size(runData.R2full_te{timeSeg},1))
            meanCur = nanmean(runData.R2full_te{timeSeg}(k,:),2);
            if meanCur >  passThreshold
                passROIs = passROIs + 1;
                
%                 et = unique(runData.eventsTypes);
%                 for z = 1:length(et)
%                     contR = nanmean(runData.R2p_test{timeSeg}(k, z, :));
%                     if ~isnan(contR)
%                         resultsByCont.(tableResults.Type{i}).(et{z}).v = resultsByCont.(tableResults.Type{i}).(et{z}).v + nanmean(runData.R2p_test{timeSeg}(k, z, :));
%                         resultsByCont.(tableResults.Type{i}).(et{z}).c = resultsByCont.(tableResults.Type{i}).(et{z}).c + 1; 
%                     end
%                 end
            end
            
            resultsByR2.(tableResults.Type{i}).r(end+1) = meanCur;            
        end
        
        resultsByAnimal.(sprintf('Animal%s', tableResults.Animal{i})).(tableResults.Neuron{i}).(tableResults.Type{i}).sum(end+1) = passROIs / totalROIs;
    end
    
    classesType = unique(tableResults.Type);
    
    for cl = 1:length(classesType)
        f1 = figure;
        hold on;

        animalsFields = fieldnames(resultsByAnimal);
        namesLabels = {};
        summaryArray = [];
        for i = 1:length(animalsFields)
            neuronNames = fieldnames(resultsByAnimal.(animalsFields{i}));
            for k = 1:length(neuronNames)
                if isfield(resultsByAnimal.((animalsFields{i})).(neuronNames{k}), classesType{cl})
                    namesLabels(end+1) = {sprintf('%s-%s', animalsFields{i}, neuronNames{k})};
                    summaryArray(end+1) = mean(resultsByAnimal.((animalsFields{i})).(neuronNames{k}).(classesType{cl}).sum);
                end
            end
        end
        
        h1 = histogram(summaryArray, 'DisplayName', classesType{cl}, 'BinWidth', 0.1);
        h1.Normalization = 'probability';
        xlim([0,1]);
        
%         bar(categorical(namesLabels), summaryArray);
        title({'Precentage Of Significant ROI"s Full Model Per Neuron', classesType{cl}});
        xlabel('ROI"s precentage pass 0.15 significant');
        ylabel('Neurons probability'); 
        mysave(f1, [outputpath, '\SignificantPrecentagePerNeuron_', classesType{cl}]);
    end

    f2 = figure;
    hold on;
    for cl = 1:length(classesType)
        h1 = histogram(resultsByR2.(classesType{cl}).r, 'DisplayName', classesType{cl});
        h1.Normalization = 'cdf';
        hold on;
    end
    xlabel('R2');
    xlim([0,1]);
    legend('show');
    
    mysave(f2, [outputpath, '\ProbabilityR2']);
    
    indexIncludeR2Big =  tableResultsAll.includeForR2 == 1 & tableResultsAll.classification == 1;
    indexIncludeMantelBig =  tableResultsAll.includeForMantel == 1 & tableResultsAll.classification == 1;
    
    r2meanBig = mean(tableResultsAll.R2All(indexIncludeR2Big));
    r2stdBig = std(tableResultsAll.R2All(indexIncludeR2Big));
  
    mantelmeanBig = mean(tableResultsAll.SubMantelValue(indexIncludeMantelBig));
    mantelstdBig = std(tableResultsAll.SubMantelValue(indexIncludeMantelBig));
  
    
    indexIncludeR2Small =  tableResultsAll.includeForR2 == 1 & tableResultsAll.classification == 2;
    indexIncludeMantelSmall =  tableResultsAll.includeForMantel == 1 & tableResultsAll.classification == 2;
    
    r2meanSmall = mean(tableResultsAll.R2All(indexIncludeR2Small));
    r2stdSmall = std(tableResultsAll.R2All(indexIncludeR2Small));
    
    mantelmeanSmall = mean(tableResultsAll.SubMantelValue(indexIncludeMantelSmall));
    mantelstdSmall = std(tableResultsAll.SubMantelValue(indexIncludeMantelSmall));
  
    f = figure;
    hold on;
    boxchart([tableResultsAll.classification(indexIncludeMantelBig);tableResultsAll.classification(indexIncludeMantelSmall)],...
        [tableResultsAll.SubMantelValue(indexIncludeMantelBig); tableResultsAll.SubMantelValue(indexIncludeMantelSmall)]);
    title('Mantel Value Small vs Big');
    plot(1:2, [mantelmeanBig, mantelmeanSmall],'--*');
    [p,h] = ranksum(tableResultsAll.SubMantelValue(indexIncludeMantelBig), tableResultsAll.SubMantelValue(indexIncludeMantelSmall));
    fileID = fopen([outputpathAll,'\MantelSummary.txt'],'w');
    fprintf(fileID,'Big mean: %f, stf: %f \n',mantelmeanBig, mantelstdBig);
    fprintf(fileID,'Small mean: %f, stf: %f \n',mantelmeanSmall, mantelstdSmall);
    fprintf(fileID,'Ttest2 h=%f, p=%f \n',h,p);
    fclose(fileID);
    
    mysave(f, [outputpathAll,'\MantelSummaryPlot']);
   
    f = figure;
    hold on;
    boxchart([tableResultsAll.classification(indexIncludeR2Big);tableResultsAll.classification(indexIncludeR2Small)],...
        [tableResultsAll.R2All(indexIncludeR2Big); tableResultsAll.R2All(indexIncludeR2Small)]);
    title('R2 Small vs Big');
    plot(1:2, [r2meanBig, r2meanSmall],'--*');
    [p,h] = ranksum(tableResultsAll.R2All(indexIncludeR2Big), tableResultsAll.R2All(indexIncludeR2Small));
    fileID = fopen([outputpathAll,'\R2Summary.txt'],'w');
    fprintf(fileID,'Big mean: %f, stf: %f \n',r2meanBig, r2stdBig);
    fprintf(fileID,'Small mean: %f, stf: %f \n',r2meanSmall, r2stdSmall);
    fprintf(fileID,'ranksum h=%f, p=%f \n',h,p);
    fclose(fileID);
    
    mysave(f, [outputpathAll,'\R2SummaryPlot']);

    
    slopeCalc = zeros(size(tableResultsAll,1),1);
    for i = 1:size(tableResultsAll,1)
        x = [];
        y = [];
        if ~isempty(tableResultsAll.GLMVectorLocation{i})
            filesL = dir(fullfile(tableResultsAll.GLMVectorLocation{i}, 'AVSDForROI_*_Depth2_*.csv'));
            
            for k = 1:length(filesL)
                data1 = readtable(fullfile(filesL(k).folder, filesL(k).name));
                x(end+1:end+size(data1, 1)) = data1.Var1;
                y(end+1:end+size(data1, 1)) = data1.Var2;
            end
            
            mdAll1 = fitglm(x ./ max(x), y);
            b1_cluster1 = mdAll1.Coefficients.Estimate(2);

            slopeCalc(i) = b1_cluster1;
        else
            slopeCalc(i) = nan;
        end
    end
    
    f = figure;
    hold on;
    boxchart([tableResultsAll.classification(indexIncludeR2Big);tableResultsAll.classification(indexIncludeR2Small)],...
        [slopeCalc(indexIncludeR2Big); slopeCalc(indexIncludeR2Small)]);
    title('Slop Value Small vs Big');
    
    plot(1:2, [mean(slopeCalc(indexIncludeR2Big)), mean(slopeCalc(indexIncludeR2Small))],'--*');
    [p,h] = ranksum(slopeCalc(indexIncludeR2Big), slopeCalc(indexIncludeR2Small));
    fileID = fopen([outputpathAll,'\SlopSummary.txt'],'w');
    fprintf(fileID,'Big mean: %f, stf: %f \n',mean(slopeCalc(indexIncludeR2Big)), std(slopeCalc(indexIncludeR2Big)));
    fprintf(fileID,'Small mean: %f, stf: %f \n',mean(slopeCalc(indexIncludeR2Small)), std(slopeCalc(indexIncludeR2Small)));
    fprintf(fileID,'ranksum h=%f, p=%f \n',h,p);
    fclose(fileID);
    
    mysave(f, [outputpathAll,'\SlopSummaryPlot']);
    
end