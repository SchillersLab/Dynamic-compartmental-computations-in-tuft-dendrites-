function mainRunnerNeuronTreeAndActivityAnalysis_HS
%     fileToRun = {'all', 'c0', 'c1', 'c2', 'c3', 'c4'};
    fileToRun = {'all', 'c0', 'c1', 'c2', 'c3'};
    
    %     Can be Euclidean OR ShortestPath OR HyperbolicDist_L OR HyperbolicDist_P  = ( Between 2 roi according to the tree path )
    roiTreeDistanceFunction = {'HyperbolicDist_P' , 'ShortestPath'};
    structureType = {'HS', 'SP'};
    
    for k = 1:length(structureType)
        for i = 1:length(fileToRun)
            HS_RunnerForCluster(fileToRun{i}, structureType{k}, roiTreeDistanceFunction{k}, 'ByH');
        end
    end
end

function HS_RunnerForCluster(clusterName, structureType, roiTreeDistanceFunction, clusterType)
    MainFolder = 'C:\Users\Jackie\Dropbox (Technion Dropbox)\Yara\Layer 5_Analysis\';
    AnimalName = 'SM04';
    DateAnimal = '08_18_19_tuft_Final_Version';
    swcFile = '08.18.19_Tuft_Final_Version_07.15.20\swcFiles_neuron2\neuron_2.swc';
    neuronNumberName = 'N2';
    RunnerDate = '2-11-20';
    RunnerNumber = 'Run1';
    behaveType = 'no_behave';
   
    neuronTreePathSWC = sprintf('%s\\Shay\\%s\\%s\\%s', MainFolder, AnimalName, DateAnimal, swcFile);
    
    neuronActivityMatrix = sprintf('%s\\Shay\\%s\\%s\\Analysis\\%s\\Structural_VS_Functional\\%s\\%s\\HS_create\\%s\\%s\\%s\\matlab_matrix_0_normal_500_100000.mat', MainFolder, AnimalName, DateAnimal, neuronNumberName, RunnerDate, RunnerNumber, behaveType, clusterType, clusterName);
    
    hyperbolicDistMatrixLocation = sprintf('%s\\Shay\\%s\\%s\\Analysis\\%s\\Structural_VS_Functional\\%s\\%s\\HS_create\\StructureTreeHS\\matlab_matrixbernoulli_100_3000.mat', MainFolder, AnimalName, DateAnimal, neuronNumberName, RunnerDate, RunnerNumber); 
    
    behaveFileTreadMillPath = '';
    behaveFrameRate = 100;
    
    outputpath = sprintf('%s\\Shay\\%s\\%s\\Analysis\\%s\\Structural_VS_Functional\\%s\\%s\\%s\\HSActivity\\%s\\%s\\', MainFolder, AnimalName, DateAnimal, neuronNumberName, RunnerDate, RunnerNumber, behaveType, structureType, clusterName);
    outputpath = char(outputpath);
    mkdir(outputpath);     
    
    excludeRoi = [];
    apical_roi = [];
    
    doComboForCloseRoi = false;
    
    firstDepthCompare = 1;
    secDepthCompare = 2;
    
%     load Tree Data
    [gRoi, ~, selectedROITable] = loadSwcFile(neuronTreePathSWC, outputpath, doComboForCloseRoi);
    
    for i = 1:length(excludeRoi)
        ex_results = contains(selectedROITable.Name, sprintf('roi%05d', excludeRoi(i)));
        
        if sum(ex_results) == 1
            selectedROITable(ex_results, :) = [];
        end
    end
    
    selectedROI = selectedROITable.Name;
    
    index_apical = zeros(1, length(apical_roi));   
    for i = 1:length(apical_roi)
        ex_results = find(contains(selectedROITable.Name, sprintf('roi%05d', apical_roi(i))));
        
        if ~isempty(ex_results)
            index_apical(i) = ex_results;
        end
    end
 
    
%     Behave TreadMillData
    if ~isempty(behaveFileTreadMillPath)
        [speedBehave, accelBehave, speedActivity, accelActivity] = treadmilBehave(behaveFileTreadMillPath, behaveFrameRate);
        plotBaseTreadMillActivity(speedBehave, accelBehave, roiActivity, outputpath);
    end
      
   if ~strcmp(hyperbolicDistMatrixLocation, "") 
        load(hyperbolicDistMatrixLocation);
   end
    
%     Calc Distance Matrix for ROI in Tree
    switch(roiTreeDistanceFunction)
        case 'Euclidean'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_Euclidean(gRoi, selectedROI, outputpath); 
        case 'ShortestPath'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_ShortestPath(gRoi, selectedROITable, outputpath);
        case 'HyperbolicDist_L'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_Hyperbolic(gRoi, selectedROITable, outputpath, loranzDistMat);
        case 'HyperbolicDist_P'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_Hyperbolic(gRoi, selectedROITable, outputpath, poincareDistMat);
    end
    
%     Calc branching
    selectedROISplitDepth1 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth1 = getSelectedROISplitBranchID(gRoi, firstDepthCompare, selectedROISplitDepth1, selectedROI);   
  
    selectedROISplitDepth3 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth3 = getSelectedROISplitBranchID(gRoi, secDepthCompare, selectedROISplitDepth3, selectedROI);   
        
    
%     Calc Activity Events Window
    snapnow
    close all;
   
   load(neuronActivityMatrix);
   
   plotTreeAndActivityDendogram(outputpath, 'HS', roiActivityDistanceMatrix, selectedROI, roiLinkage,  roiSortedByCluster, true, true);
        
   plotTreeAndActivityDendogram(outputpath, 'HS', roiActivityDistanceMatrix, selectedROI, roiLinkage,  roiSortedByCluster, true, false);
   
   plotResultesByClusterType(roiActivityDistanceMatrix, selectedROI, roiSortedByCluster, outputpath,...
       gRoi, roiTreeDistanceMatrix, selectedROISplitDepth3, selectedROISplitDepth1, selectedROITable, index_apical, clusterName, clusterType);

   snapnow
   fclose('all');
   close all;


end
 
function index_presentaion = plotResultesByClusterType(roiActivityDistanceMatrix, selectedROI, roiSortedByCluster, outputpathCurr, ...
    gRoi, roiTreeDistanceMatrix, selectedROISplitDepth3, selectedROISplitDepth1, selectedROITable, index_apical, clusterName, clusterType)
       
        figDist = figure;
        hold on;
        title({'ROI Activity Distance'});
        xticks(1:length(selectedROI));
        yticks(1:length(selectedROI));
        m = imagesc(roiActivityDistanceMatrix(roiSortedByCluster, roiSortedByCluster));
        colorbar
        cmap = jet();
        
        cmap = flipud(cmap);
        colormap(cmap);
        
        set(m,'AlphaData',~isnan(roiActivityDistanceMatrix(roiSortedByCluster, roiSortedByCluster)))
                
        for index_roi = 1:length(selectedROI)
            labelsNames(index_roi) = {sprintf('roi%d', sscanf(selectedROI{index_roi}, 'roi%d'))};
        end
       
        
        xticklabels(labelsNames(roiSortedByCluster));
        xtickangle(90);
        yticklabels(labelsNames(roiSortedByCluster));
        picNameFile = [outputpathCurr, '\DistMatrixActivity_HS'];
        mysave(figDist, picNameFile);  
        
        y = squareform(roiActivityDistanceMatrix);
        l = linkage(y, 'single');

        figDendrogram = figure;
        leafOrder = optimalleaforder(l,y);
        dendrogram(l, 'Labels', selectedROI, 'Reorder', leafOrder);
        xtickangle(90);
        
        mysave(figDendrogram, [outputpathCurr, '\DendrogramROIActivity']);
    
        plotROIDistMatrixTreeVSActivity([],gRoi, [outputpathCurr],selectedROISplitDepth1, selectedROISplitDepth1, roiTreeDistanceMatrix, roiActivityDistanceMatrix, true, 'HS', clusterName, selectedROI,index_apical,  'Distance', clusterType);

        plotROIDistMatrixTreeVSActivity([],gRoi, [outputpathCurr],selectedROISplitDepth1, selectedROISplitDepth3, roiTreeDistanceMatrix, roiActivityDistanceMatrix, false, 'HS', clusterName, selectedROI, index_apical, 'Distance', clusterType);
        
        plotRoiDistMatrixTreeVsActivityForDepthCompare(gRoi, [outputpathCurr],  selectedROITable, roiTreeDistanceMatrix, roiActivityDistanceMatrix, 'HS', clusterName);
      
end