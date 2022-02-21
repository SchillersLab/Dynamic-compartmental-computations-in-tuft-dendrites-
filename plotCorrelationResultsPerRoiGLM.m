function plotCorrelationResultsPerRoiGLM(roiActivityDistanceMatrixByH, selectedROI, selectedROISplitDepth1, colorMatrix1, outputpath)
    indexList = 1:length(selectedROI);
    meanCorr = zeros(1, length(selectedROI));
    for i = indexList
        corrV = roiActivityDistanceMatrixByH(i, indexList ~= i & selectedROISplitDepth1' == selectedROISplitDepth1(i));
        meanCorr(i) = mean(corrV);
    end
    
    classesU = unique(selectedROISplitDepth1);
    
    f = figure;
    hold on;
    title('Mean Pearson Values Per ROI, within tree');
    xlabel('Pearson');
    for i = 1:length(classesU)
        colorCur = colorMatrix1(find(selectedROISplitDepth1 == classesU(i), 1, 'first'), :);
        histogram(meanCorr(selectedROISplitDepth1 == classesU(i)),'FaceColor',colorCur, 'BinWidth', 0.1, 'FaceAlpha', 0.3);
    end
    
    mysave(f, [outputpath, '\MeanPearsonWithinTree']);
    t = table(meanCorr', selectedROI, 'VariableNames', {'RoiName', 'Mean_Correlation_Within_Tree'});
    writetable(t, [outputpath, '\MeanPearsonWithinTreeSummaryT.csv']);
end