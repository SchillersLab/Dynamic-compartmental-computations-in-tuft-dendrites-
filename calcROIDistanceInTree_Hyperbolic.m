function [roiTreeDistanceMatrix, roiSortedByCluster, l] = calcROIDistanceInTree_Hyperbolic(gRoi, selectedROI, outputpath, loranzORPDistMat, selectedROISplitDepth1, reverseHeatMap)
   tickLabels = [];
   loranzLocation = zeros(1, length(selectedROI.ID));
   
   for index = 1:length(selectedROI.ID)
        tickLabels{index} = selectedROI.Name{index};
        loranzLocation(index) = find(gRoi.Nodes.ID == selectedROI.ID(index));
   end
      
   roiTreeDistanceMatrix = loranzORPDistMat(loranzLocation, loranzLocation); 
   
    y = squareform(roiTreeDistanceMatrix);
    l = linkage(y, 'single');
         
    figDendrogram = figure;
    
    leafOrder = optimalleaforder(l,y);
    
    dendrogram(l, 'Labels', tickLabels, 'reorder', leafOrder);
    xtickangle(90);
    title('Tree Structure Dendrogram');
    mysave(figDendrogram, [outputpath, '\DendrogramROIHSDist']);
    
    roiSortedByCluster = leafOrder;    

    if selectedROISplitDepth1(roiSortedByCluster(1)) > min(selectedROISplitDepth1)
        roiSortedByCluster = roiSortedByCluster(end:-1:1);
    end
   
    if  reverseHeatMap
        roiSortedByCluster = roiSortedByCluster(end:-1:1);
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
    
    for index_roi = 1:length(tickLabels)
        labelsNames(index_roi) = {sprintf('roi%d', sscanf(tickLabels{index_roi}, 'roi%d'))};
    end
 
    
    xticklabels(labelsNames(roiSortedByCluster));
    xtickangle(90);
    yticklabels(labelsNames(roiSortedByCluster));
    
    mysave(figDist, [outputpath, '\DistMatrixROIStructure']);
end