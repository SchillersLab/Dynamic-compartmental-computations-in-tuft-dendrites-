function mainRunnerMentalTest(globalParameters)
    outputpath = globalParameters.outputpath;
    neuronTreePathSWC = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName, globalParameters.DateAnimal, globalParameters.swcFile);
    
    [gRoi, rootNodeID, selectedROITable] = loadSwcFile(neuronTreePathSWC, outputpath, false);
    
    for i = 1:length(globalParameters.excludeRoi)
        ex_results = contains(selectedROITable.Name, sprintf('roi%05d', globalParameters.excludeRoi(i)));
        
        if sum(ex_results) == 1
            selectedROITable(ex_results, :) = [];
        end
    end
    
    mkdir(fullfile(outputpath, 'mental'));
    
    selectedROI = selectedROITable.Name;

    selectedROISplitDepth1 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth1 = getSelectedROISplitBranchID(gRoi, globalParameters.firstDepthCompare, selectedROISplitDepth1, selectedROI);   
 
    selectedROISplitDepth2 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth2 = getSelectedROISplitBranchID(gRoi, globalParameters.secDepthCompare, selectedROISplitDepth2, selectedROI);   
 
    
    selectedROISplitDepth3 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth3 = getSelectedROISplitBranchID(gRoi, globalParameters.thDepthCompare, selectedROISplitDepth3, selectedROI);   

    
    switch(globalParameters.roiTreeDistanceFunction)
        case 'Euclidean'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_Euclidean(gRoi, selectedROI, outputpath, selectedROISplitDepth1); 
        case 'ShortestPath'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_ShortestPath(gRoi, selectedROITable, outputpath, selectedROISplitDepth1);
        case 'HyperbolicDist_L'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_Hyperbolic(gRoi, selectedROITable, outputpath, loranzDistMat, selectedROISplitDepth1);
        case 'HyperbolicDist_P'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_Hyperbolic(gRoi, selectedROITable, outputpath, poincareDistMat, selectedROISplitDepth1);
    end
        
    
    roiTreeDistanceMatrix = roiTreeDistanceMatrix + (diag(ones(1, size(roiTreeDistanceMatrix, 1))));
    
    for i_cluster = -1:globalParameters.clusterCount
        if i_cluster == 0
            roiActivityPeakSize = 'All';          
        elseif i_cluster == -1
            roiActivityPeakSize = 'All_ExcludeBigEvents';      
        else
            roiActivityPeakSize = ['cluster', num2str(i_cluster)];
        end
        
        outputpathCurr = [outputpath, '\', roiActivityPeakSize];
        load([outputpathCurr, '\ActivityMatrixForHS.mat'], 'roiActivityDistanceMatrixByH');
      
%         roiTreeDistanceMatrix = max(roiTreeDistanceMatrix, [], 'all') - roiTreeDistanceMatrix;
%         roiTreeDistanceMatrix = 1 - (roiTreeDistanceMatrix ./ max(roiTreeDistanceMatrix, [], 'all'));
%                  
%         first comapre all the matrix values 
        
        [pval_all(1, 1), pval_all(1, 2)] = bramila_mantel(roiActivityDistanceMatrixByH,roiTreeDistanceMatrix,1000,'pearson');
        
        NamesTree(1) = {'All'};
        
       
%         sec compare 1depth subtree
        index_p = 2;
        classes1 = unique(selectedROISplitDepth1);
        for i = 1:length(classes1)
            location = selectedROISplitDepth1 == classes1(i);
            if sum(location) <= 2
                continue;
            end
           
            [pval_all(index_p, 1), pval_all(index_p, 2)] = bramila_mantel(roiActivityDistanceMatrixByH(location, location),...
                roiTreeDistanceMatrix(location, location),1000,'pearson');  
            
            if classes1(i) == -1
                NamesTree(end + 1) = {'Depth1_NotInDepth'};   
            else
                NamesTree(end + 1) = {sprintf('Depth1_%s', gRoi.Nodes.Name{classes1(i)})};       
            end
            
            index_p = index_p + 1;
        end
%         th compare 2depth subtree
        
        classes2 = unique(selectedROISplitDepth2);
        for i = 1:length(classes2)
            location = selectedROISplitDepth2 == classes2(i);
            
            if sum(location) <= 2
                continue;
            end
            
           [pval_all(index_p, 1), pval_all(index_p, 2)] = bramila_mantel(roiActivityDistanceMatrixByH(location, location),...
                roiTreeDistanceMatrix(location, location),1000,'pearson');  
            
            if classes2(i) == -1
                NamesTree(end + 1) = {'Depth2_NotInDepth'};   
            else
                NamesTree(end + 1) = {sprintf('Depth2_%s', gRoi.Nodes.Name{classes2(i)})};       
            end
            
            index_p = index_p + 1;
        end
        
%         4 compare 3depth subtree
        
        classes3 = unique(selectedROISplitDepth3);
        for i = 1:length(classes3)
            location = selectedROISplitDepth3 == classes3(i);
            
            if sum(location) <= 2
                continue;
            end
           
            [pval_all(index_p, 1), pval_all(index_p, 2)] = bramila_mantel(roiActivityDistanceMatrixByH(location, location),...
                roiTreeDistanceMatrix(location, location),1000,'pearson');
            if classes3(i) == -1
                NamesTree(end + 1) = {'Depth3_NotInDepth'};   
            else
                NamesTree(end + 1) = {sprintf('Depth3_%s', gRoi.Nodes.Name{classes3(i)})};       
            end
            
            index_p = index_p + 1;
        end
 
        
%         save Results 
        tR = table(pval_all);
        tR.Properties.VariableNames = {'pval'};
        tR.Properties.RowNames = NamesTree;
        
        mkdir(fullfile(outputpath, 'mental', roiActivityPeakSize));
        
        writetable(tR, fullfile(outputpath, 'mental', roiActivityPeakSize, 'ResultsTable.csv'));
        
        pval_all = [];
        NamesTree = {};
%         plot Results
%         fig = uifigure;
%         uitable(fig,'Data',tR);
%         
%         mysave(fig, fullfile(outputpath, 'mental', roiActivityPeakSize, 'ResultsTableFig'));
    end
    
end