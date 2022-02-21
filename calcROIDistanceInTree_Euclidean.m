function  [roiTreeDistanceMatrix, roiSortedByCluster, l] = calcROIDistanceInTree_Euclidean(gRoi, selectedROI, outputpath, selectedROISplitDepth1, reverseHeatMap)
    roiTreeDistanceMatrix = zeros(length(selectedROI), length(selectedROI));
    for index = 1:length(selectedROI)
        for secIndex = 1:length(selectedROI)
            fNode = gRoi.Nodes(findnode(gRoi,selectedROI{index}), :);
            sNode = gRoi.Nodes(findnode(gRoi,selectedROI{secIndex}), :);
  
            roiTreeDistanceMatrix(index, secIndex) = norm([fNode.X(1), fNode.Y(1), fNode.Z(1)] - [sNode.X(1), sNode.Y(1), sNode.Z(1)]); 
        end
    end
    
    clusterNumber = sum(gRoi.Nodes.Depth(: ,1) == 2) * 3;
    
    y = squareform(roiTreeDistanceMatrix);
    l = linkage(y, 'single');
    
    leafOrder = optimalleaforder(l,y);
          
    figDendrogram = figure;
    dendrogram(l,size(roiTreeDistanceMatrix, 1), 'Labels', selectedROI, 'Reorder', leafOrder);
    xtickangle(90);
    
    title('Tree Structure Dendrogram');
    mysave(figDendrogram, [outputpath, '\DendrogramROIEuclideanDist']);
    
    roiSortedByCluster = leafOrder;
    
    if selectedROISplitDepth1(roiSortedByCluster(1)) > min(selectedROISplitDepth1)
        roiSortedByCluster = roiSortedByCluster(end:-1:1);
    end
    
       
    if  reverseHeatMap
        roiSortedByCluster = roiSortedByCluster(end:-1:1);
    end
    
    for index_roi = 1:length(selectedROI)
        labelsNames(index_roi) = {sprintf('roi%d', sscanf(selectedROI{index_roi}, 'roi%d'))};
    end
   
    figDist = figure;
    hold on;
    title({'ROI Structure Distance'});
    xticks(1:length(selectedROI));
    yticks(1:length(selectedROI));
    imagesc(roiTreeDistanceMatrix(roiSortedByCluster, roiSortedByCluster));
    colorbar
    colormap(jet);
    colormap(flipud(jet));
    
    xticklabels(labelsNames(roiSortedByCluster));
    xtickangle(90);
    yticklabels(labelsNames(roiSortedByCluster));
    
    mysave(figDist, [outputpath, '\DistMatrixROIStructure']);    
end