function [roiTreeDistanceMatrix, roiSortedByCluster, l] = calcROIDistanceInTree_Branch(gRoi, selectedROI, outputpath, selectedROISplitDepth1, reverseHeatMap)
    roiTreeDistanceMatrix = zeros(length(selectedROI.ID), length(selectedROI.ID));
    tickLabels = [];
    
   for index = 1:length(selectedROI.ID)
        tickLabels{index} = selectedROI.Name{index};
        for secIndex = 1:length(selectedROI.ID)
           
           [p, ~] = shortestpath(gRoi,selectedROI.ID(index),selectedROI.ID(secIndex), 'Method', 'positive');
           countCost = 0;
           for p_i = 1:length(p)-1
               if gRoi.Nodes.Depth(p(p_i), 1) ~= gRoi.Nodes.Depth(p(p_i + 1), 1)
                    countCost = countCost + 1;
               end
           end
           
           if countCost == 0
               roiTreeDistanceMatrix(index, secIndex) = 0;
           else
               roiTreeDistanceMatrix(index, secIndex) = (countCost - 1);
           end
           
        end
    end
      
    y = squareform(roiTreeDistanceMatrix);
    l = linkage(y, 'single');
         
    figDendrogram = figure;
    
    leafOrder = optimalleaforder(l,y);
    
    dendrogram(l, size(roiTreeDistanceMatrix, 1), 'Labels', tickLabels, 'reorder', leafOrder);
    xtickangle(90);
    title('Tree Structure Dendrogram');
    mysave(figDendrogram, [outputpath, '\DendrogramROITreeDist']);
    
    roiSortedByCluster = leafOrder;    

    if selectedROISplitDepth1(roiSortedByCluster(1)) > min(selectedROISplitDepth1)
        roiSortedByCluster = roiSortedByCluster(end:-1:1);
    end
    
    if  reverseHeatMap
        roiSortedByCluster = roiSortedByCluster(end:-1:1);
    end
    
    for index_roi = 1:length(tickLabels)
        labelsNames(index_roi) = {sprintf('roi%d', sscanf(tickLabels{index_roi}, 'roi%d'))};
    end
 
    
    figDist = figure;
    hold on;
    title({'ROI Structure Distance'});
    xticks(1:length(selectedROI.ID));
    yticks(1:length(selectedROI.ID));
    imagesc(roiTreeDistanceMatrix(roiSortedByCluster, roiSortedByCluster));
    colorbar
    colormap(jet);
    colormap(flipud(jet));
    xticklabels(labelsNames(roiSortedByCluster));
    xtickangle(90);
    yticklabels(labelsNames(roiSortedByCluster));
    
    mysave(figDist, [outputpath, '\DistMatrixROIStructure']);
end