function claculateSPMatrix(outputPath, swcFileLocation, roisTobeExcluded)
    [gRoi, rootNodeID, selectedROITable] = loadSwcFile(swcFileLocation, outputPath, false);
    
    for i = 1:length(roisTobeExcluded)
        ex_results = contains(selectedROITable.Name, sprintf('roi%05d', globalParameters.excludeRoi(i)));
        
        if sum(ex_results) == 1
            selectedROITable(ex_results, :) = [];
        end
    end
    
    selectedROI = selectedROITable.Name;
    
    selectedROISplitDepth1 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth1 = getSelectedROISplitBranchID(gRoi, 1, selectedROISplitDepth1, selectedROI, rootNodeID);   
  
    
    [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_ShortestPath(gRoi,...
        selectedROITable, outputPath, selectedROISplitDepth1, false);
     
    mysave([outputPath, '\ResultsSPCalc.mat'], 'gRoi', 'selectedROITable', 'selectedROISplitDepth1', 'SPMatrix', 'roiSortedByCluster', 'roiLinkage');
end