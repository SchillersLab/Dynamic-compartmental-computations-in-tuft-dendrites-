function mainRunnerNeuronTreeAndActivityAnalysis
    neuronTreePathSWC = "E:\Dropbox (Technion Dropbox)\Yara\Analysis\Yara's Data For Shay\SM04\08.18.19_Tuft_new_05.04.20\swcFiles\neuron_1.swc";
    
    activityByCSV = false;
    neuronActiityPathCSV = 'D:\Shay\work\N\SM04\2019-08-18_SM04_handreach_aligned_Ch1_dffByRoi.csv';
    neuronActivityPathTPA = "E:\Dropbox (Technion Dropbox)\Yara\Analysis\Yara's Data For Shay\SM04\08.18.19_Tuft_new_05.04.20\18_08_19";
    
    outputpath = "E:\Dropbox (Technion Dropbox)\Test\BYMLSPIKE";
    outputpath = char(outputpath);
    
%     Can be Euclidean OR ShortestPath = ( Between 2 roi according to the tree path )
    roiTreeDistanceFunction = 'ShortestPath';
    
%     Can be Euclidean OR PeakDistance OR WindowEventFULLPearson OR
%     WindoEventToPeakPearson OR PeaksPearson
    roiActivityDistanceFunction = 'WindoEventToPeakPearson';
        
    threshold_std = 3;
    PeakWidth = 3;
    clusterPrecentage = 4;
    histBinWidth = 0.01;
    
    roi_count = 20;
    aV = ones(1, roi_count)*0.4;
    aV(6) = 0.1;
    aV(7) = 0.1;
    aV(8) = 0.05;
    aV(9) = 0.1;
    aV(11) = 0.1;
    aV(12) = 0.1;
    aV(14) = 0.2;
    aV(15) = 0.3;
    aV(16) = 0.1;
    aV(17) = 0.1;
    aV(18) = 0.1;
    
%     Trial Number To plot with Tree
    trialNumber = [1, 2] ;
    samplingRate = 20;
    time_sec_of_trial = 12;
    
    
    firstDepthCompare = 2;
    secDepthCompare = 3;
    
%     load Tree Data
    [gRoi, ~, selectedROITable] = loadSwcFile(neuronTreePathSWC, outputpath);
    
    selectedROI = selectedROITable.Name;
    
    selectedROISplitDepth1 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth1 = getSelectedROISplitBranchID(gRoi, firstDepthCompare, selectedROISplitDepth1, selectedROI);   
    
    if (activityByCSV)
        %     load roi activity file
        [roiActivity, roiActivityNames] = loadActivityFile(neuronActiityPathCSV, selectedROI);
    else
        [roiActivity, roiActivityNames, tr_frame_count] = loadActivityFileFromTPA(neuronActivityPathTPA, selectedROI, outputpath);
    end
        
    
%     Calc Distance Matrix for ROI in Tree
    switch(roiTreeDistanceFunction)
        case 'Euclidean'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_Euclidean(gRoi, selectedROI, outputpath); 
        case 'ShortestPath'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_ShortestPath(gRoi, selectedROITable, outputpath);
    end
    
%     Calc Activity Events Window
%     [windowFULLOrig, windowToPeakOrig, locationPeaksOrig, activityClusterValue, clusterCount] = calcActivityEventsWindowsAndPeaks(roiActivity, outputpath, threshold_std, PeakWidth, clusterPrecentage, histBinWidth, samplingRate);
    
    [windowFULLOrig, windowToPeakOrig, locationPeaksOrig, activityClusterValue, clusterCount] = calcActivityEventsWindowsAndPeaksType2(roiActivity, outputpath, threshold_std, PeakWidth, clusterPrecentage, histBinWidth, samplingRate, tr_frame_count, aV);

    
    %     Plot Tree And Trial Activity in the ROI
    
%     totalTrialTime = samplingRate * time_sec_of_trial;
%     plotTreeAndActivityForTrial(trialNumber, totalTrialTime, roiSortedByCluster, roiActivity, roiActivityNames, selectedROI, outputpath, locationPeaks, windowFULL, roiLinkage);
%     
    
    for i_cluster = 0:clusterCount
        if i_cluster == 0
            peakSelectedIndex = 1:length(locationPeaksOrig);
            roiActivityPeakSize = 'All';
        else
            peakSelectedIndex = (activityClusterValue == i_cluster);
            roiActivityPeakSize = ['cluster', num2str(i_cluster)];
        end
        
        if (length(peakSelectedIndex) < 2)
            continue;
        end
        
        outputpathCurr = [outputpath, '\', roiActivityPeakSize '\'];
           
        locationPeaks = locationPeaksOrig(peakSelectedIndex);
        windowToPeak = windowToPeakOrig(peakSelectedIndex, :);
        windowFULL = windowFULLOrig(peakSelectedIndex, :);

    %     Calc Distance Matrix for ROI in Activity
       switch(roiActivityDistanceFunction)
            case 'Euclidean'
               [roiActivityDistanceMatrix] = calcROIDistanceInActivity_Euclidean(roiActivity, roiActivityNames, selectedROI, windowFULL); 
            case 'PeakDistance'
               [roiActivityDistanceMatrix] = calcROIDistanceInActivity_PeakDistance(roiActivity, roiActivityNames, selectedROI, locationPeaks);
           case 'WindowEventFULLPearson'
               [roiActivityDistanceMatrix] = calcROIDistanceInActivity_WindowEventPearson(roiActivity, roiActivityNames, selectedROI, windowFULL);
           case 'WindoEventToPeakPearson'
               [roiActivityDistanceMatrix] = calcROIDistanceInActivity_WindowEventPearson(roiActivity, roiActivityNames, selectedROI, windowToPeak);
           case 'PeaksPearson'
               [roiActivityDistanceMatrix] = calcROIDistanceInActivity_WindowEventPearson(roiActivity, roiActivityNames, selectedROI, [locationPeaks, locationPeaks]);
       end

        figDist = figure;
        hold on;
        title({'ROI Activity Distance'});
        xticks(1:length(selectedROI));
        yticks(1:length(selectedROI));
        imagesc(roiActivityDistanceMatrix(roiSortedByCluster, roiSortedByCluster));
        colorbar
        colormap(jet);
        colormap(flipud(jet));
        
        caxis([0,1]);
        xticklabels(selectedROI(roiSortedByCluster));
        xtickangle(90);
        yticklabels(selectedROI(roiSortedByCluster));

        mysave(figDist, [outputpathCurr, '\DistMatrixActivity_', roiActivityDistanceFunction, '_eventsSize', roiActivityPeakSize]);  

    %     yAct = squareform(roiActivityDistanceMatrix);
    %     lAct = linkage(yAct, 'single');
    %           
    %     figDendrogram = figure;
    %     
    %     leafOrderAct = optimalleaforder(lAct,yAct);
    %     
    %     dendrogram(lAct, 'Labels', selectedROI, 'reorder', leafOrderAct);
    %     mysave(figDendrogram, [outputpath, '\DendrogramROIActivityDist', roiActivityDistanceFunction, '_eventsSize', roiActivityPeakSize]);
    %     

        plotROIDistMatrixTreeVSActivity(gRoi, outputpathCurr, selectedROISplitDepth1, roiTreeDistanceMatrix, roiActivityDistanceMatrix, true, roiActivityDistanceFunction, roiActivityPeakSize, selectedROI);


        selectedROISplitDepth3 = ones(length(selectedROI), 1) * -1;
        selectedROISplitDepth3 = getSelectedROISplitBranchID(gRoi, secDepthCompare, selectedROISplitDepth3, selectedROI);   

        plotROIDistMatrixTreeVSActivity(gRoi, outputpathCurr, selectedROISplitDepth3, roiTreeDistanceMatrix, roiActivityDistanceMatrix, false, roiActivityDistanceFunction, roiActivityPeakSize, selectedROI);
        
        
        plotRoiDistMatrixTreeVsActivityForDepthCompare(gRoi, outputpathCurr,  selectedROITable, roiTreeDistanceMatrix, roiActivityDistanceMatrix, roiActivityDistanceFunction, roiActivityPeakSize );
        
        valuePeaks = [];    
        locationPeaks = [];
        windowToPeak = [];
        windowFULL = [];
    end   
 end