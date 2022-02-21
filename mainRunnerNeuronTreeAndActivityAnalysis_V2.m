function mainRunnerNeuronTreeAndActivityAnalysis_V2
    neuronTreePathSWC = "C:\Users\Jackie\Dropbox (Technion Dropbox)\Yara\Layer 5_Analysis\Yara's Data For Shay\SM04\08.18.19_Tuft\08.18.19_Tuft_Final_Version_07.15.20\swcFiles_neuron1,3,4\neuron_1.swc";
    
    activityByCSV = false;
    neuronActiityPathCSV = "C:\Users\Jackie\Dropbox (Technion Dropbox)\Yara\Layer 5_Analysis\Yara's Data For Shay\SM04\08.18.19_Tuft\Nate\2019-08-18_SM04_handreach_aligned_Ch1_dffByRoi.csv";
    neuronActivityPathTPA = "C:\Users\Jackie\Dropbox (Technion Dropbox)\Yara\Layer 5_Analysis\SM04\08_18_19_tuft_Final_Version";
    
    behaveFileTreadMillPath = '';
    behaveFrameRate = 100;
    
    outputpath = "C:\Users\Jackie\Dropbox (Technion Dropbox)\Yara\Layer 5_Analysis\Yara's Data For Shay\SM04\08.18.19_Tuft\Yara\Results\Ne1_V1\";
    outputpath = char(outputpath);
    mkdir(outputpath);
   
%     Can be Euclidean OR ShortestPath = ( Between 2 roi according to the tree path )
    roiTreeDistanceFunction = 'ShortestPath';
    
%     Can be WindowEventFULLPearson OR
%     WindoEventToPeakPearson OR PeaksPearson
    roiActivityDistanceFunction = 'WindoEventToPeakPearson';
      
    excludeRoi = [9];
    
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
    
    selectedROISplitDepth1 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth1 = getSelectedROISplitBranchID(gRoi, firstDepthCompare, selectedROISplitDepth1, selectedROI);   
    
    clusterCount = 3;
    
    roi_count = length(selectedROI);
    aV = ones(1, roi_count)*0.3;      
    
    sigmaChangeValue = zeros(1, roi_count);

    DeconvFiltDur = ones(1, roi_count) * 0.4;       % smoothing filter duration in sec  
    do_roiFilter = ones(1, roi_count) * 0;
    
%     Trial Number To plot with Tree

    save([outputpath '\runParametes'],'aV', 'roi_count', 'do_roiFilter', 'DeconvFiltDur', 'sigmaChangeValue', 'excludeRoi');

    trialNumber = [1, 2] ;
    samplingRate = 20;
    time_sec_of_trial = 12;
        
    
    if (activityByCSV)
        %     load roi activity file
        [roiActivity, roiActivityNames] = loadActivityFile(neuronActiityPathCSV, selectedROI);
        tr_frame_count = [];
    else
        [roiActivity, roiActivityNames, tr_frame_count] = loadActivityFileFromTPA(neuronActivityPathTPA, selectedROI, outputpath);
    end
    
%     Behave TreadMillData
    if ~isempty(behaveFileTreadMillPath)
        [speedBehave, accelBehave, speedActivity, accelActivity] = treadmilBehave(behaveFileTreadMillPath, behaveFrameRate);
        plotBaseTreadMillActivity(speedBehave, accelBehave, roiActivity, outputpath);
    end
          
%     Calc Distance Matrix for ROI in Tree
    switch(roiTreeDistanceFunction)
        case 'Euclidean'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_Euclidean(gRoi, selectedROI, outputpath); 
        case 'ShortestPath'
           [roiTreeDistanceMatrix, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_ShortestPath(gRoi, selectedROITable, outputpath);
    end
    
    
%     Calc Activity Events Window
   
    [all_locationFull_start, all_locationFull_end, all_locationFull_pks, activityClusterValue, roiActivity_comb] = calcActivityEventsWindowsAndPeaks_V2(roiActivity, outputpath, clusterCount, samplingRate, tr_frame_count, aV, roiActivityNames, sigmaChangeValue);

    all_event_struct.end = all_locationFull_end;
    
    all_event_struct.start = all_locationFull_start;
    
    all_event_struct.pks = all_locationFull_pks;
    
    all_event_struct.cluster = activityClusterValue;
    
    %     Plot Tree And Trial Activity in the ROI
    
%     totalTrialTime = samplingRate * time_sec_of_trial;
%     plotTreeAndActivityForTrial(trialNumber, totalTrialTime, roiSortedByCluster, roiActivity, roiActivityNames, selectedROI, outputpath, locationPeaks, windowFULL, roiLinkage);
%     
    

%     Filter ROI activity if needed
%     
    sampFreq        = samplingRate;
    filtDur         = DeconvFiltDur .* sampFreq;      % filter duration in sec
    filtLenH        = ceil(filtDur/2);
    filtLen         = filtLenH*2;

    for i_activity = 1:size(roiActivity, 2)
        if do_roiFilter == 1
            filtSmooth          = hamming(filtLen(i_activity));
            filtSmooth          = filtSmooth./sum(filtSmooth);

            fig = figure;
            hold on;
            title({'ROI with filter', roiActivityNames{i_activity}});
            plot(roiActivity(:, i_activity));
            roiActivity(:, i_activity) = filtfilt(filtSmooth,1,roiActivity(:, i_activity));   
            plot(roiActivity(:, i_activity));

            mysave(fig, [outputpath '\filterResults\roiActivity_' num2str(i_activity)]); 

            close(fig);
        end
    end
    
    import mlreportgen.ppt.*
    ppt = Presentation([outputpath '\AnalysisResultsPresentation'], 'AnalysisP.potm');
    open(ppt);
    currentResultsSlide = add(ppt, 'AnalysisP');
    index_presentaion = 1;
       
    for i_cluster = -1:clusterCount
        if i_cluster == 0
            roiActivityPeakSize = 'All';
        elseif i_cluster == -1
            roiActivityPeakSize = 'All_ExcludeSmallEvents';
        else
            roiActivityPeakSize = ['cluster', num2str(i_cluster)];
        end
        
        outputpathCurr = [outputpath, roiActivityPeakSize];
           
    %     Calc Distance Matrix for ROI in Activity
       switch(roiActivityDistanceFunction)
           case 'WindowEventFULLPearson'
               [roiActivityDistanceMatrix] = calcROIDistanceInActivity_WindowEventPearson(roiActivity_comb, roiActivityNames, selectedROI, all_event_struct, 'FULL', i_cluster);
           case 'WindoEventToPeakPearson'
               [roiActivityDistanceMatrix] = calcROIDistanceInActivity_WindowEventPearson(roiActivity_comb, roiActivityNames, selectedROI, all_event_struct, 'ToPeak', i_cluster);
           case 'PeaksPearson'
               [roiActivityDistanceMatrix] = calcROIDistanceInActivity_WindowEventPearson(roiActivity_comb, roiActivityNames, selectedROI, all_event_struct, 'Peaks', i_cluster);
       end

        figDist = figure;
        hold on;
        title({'ROI Activity Distance'});
        xticks(1:length(selectedROI));
        yticks(1:length(selectedROI));
        imagesc(roiActivityDistanceMatrix(roiSortedByCluster, roiSortedByCluster));
        colorbar
        colormap(jet);
%         colormap(flipud(jet));
        
%         caxis([0,1]);
        xticklabels(selectedROI(roiSortedByCluster));
        xtickangle(90);
        yticklabels(selectedROI(roiSortedByCluster));
        picNameFile = [outputpathCurr, '\DistMatrixActivity_', roiActivityDistanceFunction, '_eventsSize', roiActivityPeakSize];
        mysave(figDist, picNameFile);  

        pictureNames = plotROIDistMatrixTreeVSActivity(gRoi, outputpathCurr, selectedROISplitDepth1, roiTreeDistanceMatrix, roiActivityDistanceMatrix, true, roiActivityDistanceFunction, roiActivityPeakSize, selectedROI);


        selectedROISplitDepth3 = ones(length(selectedROI), 1) * -1;
        selectedROISplitDepth3 = getSelectedROISplitBranchID(gRoi, secDepthCompare, selectedROISplitDepth3, selectedROI);   

        plotROIDistMatrixTreeVSActivity(gRoi, outputpathCurr, selectedROISplitDepth3, roiTreeDistanceMatrix, roiActivityDistanceMatrix, false, roiActivityDistanceFunction, roiActivityPeakSize, selectedROI);
        
        
        plotRoiDistMatrixTreeVsActivityForDepthCompare(gRoi, outputpathCurr,  selectedROITable, roiTreeDistanceMatrix, roiActivityDistanceMatrix, roiActivityDistanceFunction, roiActivityPeakSize );
        
        if i_cluster ~= -1
            replace(currentResultsSlide.Children(index_presentaion), Picture([picNameFile '.tif']));       
            replace(currentResultsSlide.Children(index_presentaion + 1), Picture([pictureNames{1} '.tif']));
            index_presentaion = index_presentaion + 2;
        end
    end  
    
    close(ppt);
 end