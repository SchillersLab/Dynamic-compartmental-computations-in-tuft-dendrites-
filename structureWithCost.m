function structureWithCost(gRoi, selectedROITable, roiActivityDistanceMatrixByH, selectedROISplitDepth1, outputpath, i_cluster, reverseHeatMap)
    costVector = 0:50:500;
    mdl = {};
    r_2 = [];
    x = [];
    indexC = 1;
    for costv = costVector
        [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_ShortestPathCost(gRoi, selectedROITable, [outputpath ,'\cost\', num2str(i_cluster), '\', num2str(costv), '\'], selectedROISplitDepth1, costv, reverseHeatMap);
        
        y = [];
        
        index_i = 1;
        for i_roi = 1:size(selectedROITable,1)
            for k_roi = (i_roi+1):size(selectedROITable,1)
                x(indexC, index_i) = roiTreeDistanceMatrix(i_roi, k_roi);
                y(index_i) = roiActivityDistanceMatrixByH(i_roi, k_roi);
                index_i = index_i + 1;
            end
        end
        
        mdl(end+1) = {fitglm(x(indexC, :), y)};
        r_2(end+1) = mdl{end}.Rsquared.Ordinary;
        
        indexC =indexC + 1;
    end
    
    calcCont = diff([r_2(1), r_2]);
    calcCont(find(calcCont < 0, 1, 'first'):end) = 0;
    
    r_2loc = find((cumsum(calcCont) ./ sum(calcCont)) >= 0.99, 1, 'first');
    
    if isempty(r_2loc)
        r_2loc = 1;
    end
    
    r2Max = r_2(r_2loc);
    
    fig = figure;
    hold on;
    plot(costVector, r_2);
    scatter(costVector(r_2loc), r2Max, 'o');
    xlabel('cost');
    ylabel('R^2');
    title('R^2 vs Cost');
    mysave(fig, [outputpath ,'\cost\', num2str(i_cluster), '\','r2Histogram']);
    
    fig = figure;
    hold on;
    scatter(x(r_2loc,:), y, 20, 'filled', 'MarkerFaceColor', [0,0,0]);
    plot(x(r_2loc,:), predict(mdl{r_2loc}, x(r_2loc,:)'), 'Color', [119, 123, 126] ./ 255, 'LineStyle', '--');
    
    pV = foldsCalculator(x(r_2loc,:), y, mdl{r_2loc}.Rsquared.Ordinary);
    
    title({sprintf('ROI Activity VS Tree Distance cost %d', costVector(r_2loc)),...
        sprintf('R^2 = %f', mdl{r_2loc}.Rsquared.Ordinary), ...
        sprintf('pvalue = %d', pV)});
    ylabel({'Calcium Events Correlation'});
    xlabel({'Dendritic distance'});
    mysave(fig, [outputpath ,'\cost\' , num2str(i_cluster), '\', 'ActivityVSDistCost']);    
end