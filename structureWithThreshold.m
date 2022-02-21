function structureWithThreshold(selectedROI,roiSortedByCluster, roiTreeDistanceMatrix, roiActivityDistanceMatrixByH, outputpath, i_cluster, tJump, roiTreeDistanceFunction, otherPath)
    minDistance = min(roiTreeDistanceMatrix, [], 'all', 'omitnan');
    min2Distance = min(roiTreeDistanceMatrix(roiTreeDistanceMatrix > 0), [], 'all', 'omitnan');
    maxDistance = max(roiTreeDistanceMatrix, [], 'all', 'omitnan');
    mkdir(outputpath);
    
    if sum(isnan(roiActivityDistanceMatrixByH), 'all') < (size(roiActivityDistanceMatrixByH(:), 1) - size(roiActivityDistanceMatrixByH, 1)) &&...
            sum((roiActivityDistanceMatrixByH < 0), 'all') < (size(roiActivityDistanceMatrixByH(:), 1) - size(roiActivityDistanceMatrixByH, 1))
       
       
       roiActivityDistanceMatrixByHNO_nan = roiActivityDistanceMatrixByH;
       roiActivityDistanceMatrixByHNO_nan(isnan(roiActivityDistanceMatrixByHNO_nan)) = 0;
       roiActivityDistanceMatrixByHNO_nan(roiActivityDistanceMatrixByHNO_nan < 0) = 0;
    
       tJM = minDistance:tJump:maxDistance;
       
       location2start = find(tJM >= min2Distance);
       
       outM(1:(location2start-1)) = 0;
       pValM(1:(location2start-1)) = 1;
       
       for i = location2start:length(tJM)
           if tJM(i) == maxDistance
              outM(i) = 0 ;
              pValM(i) = 1;
              continue;
           end
           
           dMatT = roiTreeDistanceMatrix > tJM(i);
           
           [outM(i), pValM(i)] = bramila_mantel(1 - abs(roiActivityDistanceMatrixByHNO_nan), dMatT, 5000, 'pearson');  
       end
       
       maxMantel = max(outM);
       maxI = find(outM == maxMantel, 1, 'first');
       maxPValue = pValM(maxI);
       maxThreshold = tJM(maxI);
       
%        if maxPValue > 0.05
%            maxMantel = 0;
%            maxPValue = 1;
%            maxThreshold = 0;
%        end
%        
       f = figure;
       hold on;
       scatter(tJM, outM, 'k', 'filled');
       xlabel('Threshold Value microns');
       ylabel('Mantel Value');
       mysave(f, fullfile(outputpath, sprintf('%s_cluster%d_allvaluesDist',roiTreeDistanceFunction, i_cluster)));
       
       mysave(f, fullfile(otherPath, sprintf('%s_cluster%d_allvaluesDist',roiTreeDistanceFunction, i_cluster)));
       
       
        for index_roi = 1:length(selectedROI)
            labelsNames(index_roi) = {sprintf('roi%d', sscanf(selectedROI{index_roi}, 'roi%d'))};
        end
 
        dMatTMax = roiTreeDistanceMatrix > maxThreshold;
        figDist = figure;
        hold on;
        title({sprintf('ROI Structure Distance - Threshold = %d', maxThreshold)});
        xticks(1:length(selectedROI));
        yticks(1:length(selectedROI));
        imagesc(dMatTMax(roiSortedByCluster, roiSortedByCluster));
        colorbar
        colormap(jet);
        colormap(flipud(jet));
        xticklabels(labelsNames(roiSortedByCluster));
        xtickangle(90);
        yticklabels(labelsNames(roiSortedByCluster));

        mysave(figDist, fullfile(outputpath, sprintf('%s_cluster%d_thresholdMax',roiTreeDistanceFunction, i_cluster)));

        mysave(figDist, fullfile(otherPath, sprintf('%s_cluster%d_thresholdMax',roiTreeDistanceFunction, i_cluster)));

    else
       maxMantel = 0;
       maxPValue = 1;
       maxThreshold = nan;
       outM = nan;
       pValM = nan;
    end
    
   
   save(fullfile(outputpath, sprintf('%s_cluster%d',roiTreeDistanceFunction, i_cluster)), 'maxPValue', 'maxMantel', 'maxThreshold', 'outM', 'pValM', 'min2Distance', 'maxDistance', 'tJump');
   
   save(fullfile(otherPath, sprintf('%s_cluster%d',roiTreeDistanceFunction, i_cluster)), 'maxPValue', 'maxMantel', 'maxThreshold', 'outM', 'pValM', 'min2Distance', 'maxDistance', 'tJump');
   
end