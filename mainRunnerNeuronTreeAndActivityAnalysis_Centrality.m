function mainRunnerNeuronTreeAndActivityAnalysis_Centrality(globalParameters) 
    
    neuronTreePathSWC = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName, globalParameters.DateAnimal, globalParameters.swcFile);
    
    activityByCSV = false;
    neuronActiityPathCSV = '';
    neuronActivityPathTPA = fullfile(globalParameters.TPAFolder, globalParameters.AnimalName, globalParameters.DateAnimal);
    
    outputpath = globalParameters.outputpath;
   
    behaveFileTreadMillPath = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, globalParameters.DateAnimal, globalParameters.treadmilFile);
    
    doComboForCloseRoi = false;
    
    eventsDetectionFolder = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, ...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
    globalParameters.RunnerDate,globalParameters.RunnerNumber, 'EventsDetection');
    mkdir(eventsDetectionFolder);
    
    
    centralityFolder = fullfile(globalParameters.outputpath, 'centrality');
    mkdir(centralityFolder);

    sprintf('Animal :%s, Date :%s, Neuron :%s, Behave :%s, Analysis :%s', globalParameters.AnimalName, globalParameters.DateAnimal, globalParameters.neuronNumberName, globalParameters.behaveType, globalParameters.analysisType)
    
%     load Tree Data
    sprintf('Structure Plot Results')
       
    [gRoi, rootNodeID, selectedROITable] = loadSwcFile(neuronTreePathSWC, centralityFolder, doComboForCloseRoi);
    
    for i = 1:length(globalParameters.excludeRoi)
        ex_results = contains(selectedROITable.Name, sprintf('roi%05d', globalParameters.excludeRoi(i)));
        
        if sum(ex_results) == 1
            selectedROITable(ex_results, :) = [];
        end
    end
    
    selectedROI = selectedROITable.Name;
     
    index_apical = zeros(1, length(globalParameters.apical_roi));   
    for i = 1:length(globalParameters.apical_roi)
        ex_results = find(contains(selectedROITable.Name, sprintf('roi%05d', globalParameters.apical_roi(i))));
        
        if ~isempty(ex_results)
            index_apical(i) = ex_results;
        end
    end
 
    
    roi_count = length(selectedROI);
    aV = ones(1, roi_count)*globalParameters.aVForAll;
    
    if ~isempty(globalParameters.aVFix.location)
        aV(globalParameters.aVFix.location) = globalParameters.aVFix.values;
    end
    
    sigmaChangeValue = zeros(1, roi_count);

    
%     Trial Number To plot with Tree

    save([outputpath '\runParametes'],'aV', 'roi_count', 'sigmaChangeValue', 'globalParameters');

    fid=fopen([outputpath '\Parametes.txt'],'w');
    fprintf(fid, 'hyperbolicDistMatrixLocation : %s r\n', globalParameters.hyperbolicDistMatrixLocation);    
    fprintf(fid, 'roiTreeDistanceFunction : %s r\n', globalParameters.roiTreeDistanceFunction);
    fprintf(fid, 'roiActivityDistanceFunction : %s r\n', globalParameters.roiActivityDistanceFunction);    
    fprintf(fid, 'clusterCount : %d r\n', globalParameters.clusterCount);
    fprintf(fid, 'eventWin : %d r\n', globalParameters.eventWin);
    
    for in = 1:length(globalParameters.runByEvent)
        fprintf(fid, 'event behave : %s r\n', globalParameters.runByEvent{in});
    end
    
    fprintf(fid, 'event behave lag : %d - %d r\n', globalParameters.runBehaveLag(1), globalParameters.runBehaveLag(2));
    fclose(fid);
        
    
    if (activityByCSV)
        %     load roi activity file
        [roiActivity, roiActivityNames] = loadActivityFile(neuronActiityPathCSV, selectedROI);
        tr_frame_count = [];
    else
        [roiActivity, roiActivityNames, tr_frame_count] = loadActivityFileFromTPA(neuronActivityPathTPA, selectedROI, outputpath);
    end
      
%     Behave TreadMillData
      
   if ~strcmp(globalParameters.hyperbolicDistMatrixLocation, "") 
        load(globalParameters.hyperbolicDistMatrixLocation);
   end
   
   %     Calc branching
    selectedROISplitDepth1 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth1 = getSelectedROISplitBranchID(gRoi, globalParameters.firstDepthCompare, selectedROISplitDepth1, selectedROI);   
  
    selectedROISplitDepth3 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth3 = getSelectedROISplitBranchID(gRoi, globalParameters.secDepthCompare, selectedROISplitDepth3, selectedROI);   

    % Save Graph for HS structure with colors
    classesD1 = unique(selectedROISplitDepth1);   
    classesD2 = unique(selectedROISplitDepth3);   
    
    classesD1(classesD1 == -1) = [];
    classesD2(classesD2 == -1) = [];
    
    colorMatrix1 = zeros(length(selectedROISplitDepth1), 3);
    colorMatrix2 = zeros(length(selectedROISplitDepth3), 3);
    for d_i = 1:length(selectedROISplitDepth1)
        colorMatrix1(d_i, :) = getTreeColor('within', find(classesD1 == selectedROISplitDepth1(d_i)));
        colorMatrix2(d_i, :) = getTreeColor('within', find(classesD2 == selectedROISplitDepth3(d_i)));
    end
    
%     Calc Activity Events Window
    snapnow;
    close all;

    if isfile([eventsDetectionFolder, '\roiActivity_comb.mat'])
        load([eventsDetectionFolder, '\roiActivity_comb.mat'], 'allEventsTable', 'roiActivity_comb'); 
        
        if globalParameters.reRunClusterData
        %     -----------------------------------------------------------------------------------------------------

            SpikeTrainClusterSecByH = getClusterForActivity(allEventsTable.H, globalParameters.clusterCount);
            printClusterResults(SpikeTrainClusterSecByH, globalParameters.clusterCount, mean(roiActivity_comb, 2), allEventsTable.pks, allEventsTable.start, allEventsTable.event_end, allEventsTable.H, outputpath, 'ByH')

        %   -----------------------------------------------------------------------------------------------------  

            SpikeTrainClusterSecByPrecantage = getClusterForActivity(allEventsTable.roiPrecantage, globalParameters.clusterCount);
            printClusterResults(SpikeTrainClusterSecByPrecantage, globalParameters.clusterCount, mean(roiActivity_comb, 2), allEventsTable.pks, allEventsTable.start, allEventsTable.event_end, allEventsTable.H, outputpath, 'ByP')

        %     -----------------------------------------------------------------------------------------------------

            allEventsTable.clusterByRoiPrecantage = SpikeTrainClusterSecByPrecantage';
            allEventsTable.clusterByH = SpikeTrainClusterSecByH';
        end
        
    else
        error('first run all with no events')
    end    
    
    for i_e = 1:size(allEventsTable, 1)
        tr_index = floor(allEventsTable.start(i_e) ./ tr_frame_count) + 1;
        allEventsTable.tr_index(i_e) = tr_index;
    end
    
    if globalParameters.isHandreach 
        [BehaveDataAll, NAMES, trials_label] = loadBDAFile(neuronActivityPathTPA, globalParameters.BehavioralSamplingRate, globalParameters.ImageSamplingRate, tr_frame_count, globalParameters.behavioralDelay, globalParameters.toneTime);
        
        if ~strcmp(globalParameters.excludeTrailsByEventCount.Name, 'non')
           behaveCountForCa = BehaveDataAll.(['last', globalParameters.excludeTrailsByEventCount.Name]).count(allEventsTable.tr_index);
          
           fig = figure;
           hold on;

           subplot(2, 1, 1);
           h_ByT = histogram(BehaveDataAll.(['last', globalParameters.excludeTrailsByEventCount.Name]).count); 
           title({['Event ' globalParameters.excludeTrailsByEventCount.Name], 'Histogram By Trial'});

           subplot(2, 1, 2);
           h_ByP = histogram(behaveCountForCa);
           title({['Event ' globalParameters.excludeTrailsByEventCount.Name], 'Histogram By Ca Events'});

           mysave(fig, [outputpath, '\HistogramEventsCount_' globalParameters.excludeTrailsByEventCount.Name]);  
            
            allEventsTable(behaveCountForCa < globalParameters.excludeTrailsByEventCount.countRange(1) | ...
                behaveCountForCa > globalParameters.excludeTrailsByEventCount.countRange(2),:) = [];
        end
       
       if all(~strcmp(globalParameters.runByEvent, 'non'))          
           runByEventTemp = {};
           runBehaveLagTemp = [];
           for i_run = 1:length(globalParameters.runByEvent)
                if contains(globalParameters.runByEvent(i_run), '_all')
                    runByEvent_fix = replace(globalParameters.runByEvent{i_run}, '_all', '');
                    newEvents = NAMES(contains(NAMES, runByEvent_fix)&(~contains(NAMES, ['last' runByEvent_fix]))&(~contains(NAMES, ['firstTone' runByEvent_fix])));
                    runByEventTemp((end + 1 ): (end + length(newEvents))) = newEvents;
                    runBehaveLagTemp((end + 1 ): (end + length(newEvents)), :) = [ones(length(newEvents), 1) * globalParameters.runBehaveLag(i_run, 1), ones(length(newEvents), 1) * globalParameters.runBehaveLag(i_run, 2)];
                else
                    runByEventTemp(end+1) = globalParameters.runByEvent(i_run);
                    runBehaveLagTemp(end + 1, :) = globalParameters.runBehaveLag(i_run, :);
                end
           end
           
           currentEventLoc = zeros(length(runByEventTemp),1);
          
           for in = 1:length(runByEventTemp)
              currentEventLoc(in) = find(strcmp(NAMES, runByEventTemp{in}));
              
              if ~strcmp(globalParameters.FirstEventAfter , 'non')
                   [BehaveDataAll.([NAMES{currentEventLoc(in)},'_', globalParameters.FirstEventAfter{1}]).startTiming,...
                       BehaveDataAll.([NAMES{currentEventLoc(in)},'_', globalParameters.FirstEventAfter{1}]).endTiming] = findFirstEventsAfter(BehaveDataAll, NAMES, NAMES{currentEventLoc(in)}, globalParameters.FirstEventAfter{1});
                   NAMES(end + 1) = {[NAMES{currentEventLoc(in)},'_', globalParameters.FirstEventAfter{1}]};
                  
                   currentEventLoc(in) = length(NAMES); 
                   runByEventTemp(in) =  NAMES(end);
              end
              
              if globalParameters.doBehaveAlignedPlot
                  if strcmp(globalParameters.EventTiming, 'start')
                      plotEventCaForBehaveDataHandReach(BehaveDataAll.(NAMES{currentEventLoc(in)}).startTiming, tr_frame_count, allEventsTable, globalParameters.clusterCount, outputpath, NAMES{currentEventLoc(in)}, trials_label, globalParameters.split_trialsLabel)           
                  else
                      plotEventCaForBehaveDataHandReach(BehaveDataAll.(NAMES{currentEventLoc(in)}).endTiming, tr_frame_count, allEventsTable, globalParameters.clusterCount, outputpath, NAMES{currentEventLoc(in)}, trials_label, globalParameters.split_trialsLabel)           
                  end
              end
           end
           
           eventsIndexTodelete = zeros(1, size(allEventsTable, 1));
           for i_e = 1:size(allEventsTable, 1)
               if globalParameters.split_trialsLabel ~= 0 & trials_label(allEventsTable.tr_index(i_e)) ~= globalParameters.split_trialsLabel
                    eventsIndexTodelete(i_e) = 1;
                    continue;
               end
               
               alignedLocation = zeros(1, length(runByEventTemp));
               aligned_start = zeros(1, length(runByEventTemp));
               
               checkEvent = zeros(1, length(runByEventTemp));
               for in = 1:length(runByEventTemp) 
                   if strcmp(globalParameters.EventTiming, 'start')
                       alignedLocation(in) = BehaveDataAll.(NAMES{currentEventLoc(in)}).startTiming(allEventsTable.tr_index(i_e));
                   else
                       alignedLocation(in) = BehaveDataAll.(NAMES{currentEventLoc(in)}).endTiming(allEventsTable.tr_index(i_e));
                   end
                   
                   if alignedLocation(in) ~= 0
                      alignedLocation(in) = (allEventsTable.tr_index(i_e) - 1) * tr_frame_count + alignedLocation(in);
                      aligned_start(in) = allEventsTable.start(i_e) - alignedLocation(in);
                      if (aligned_start(in) < runBehaveLagTemp(in, 1)) | (aligned_start(in) > runBehaveLagTemp(in, 2))
                         checkEvent(in) = 1; 
                      end
                   else
                       checkEvent(in) = 1;
                   end
               end               
               
                if globalParameters.do_eventsBetween
                    if ~(all(checkEvent == 0))
                       eventsIndexTodelete(i_e) = 1;
                    end
                else
                    if all(checkEvent == 1)
                        eventsIndexTodelete(i_e) = 1;
                    end
                end
           end

           
           if globalParameters.runByNegEvent
               allEventsTable(eventsIndexTodelete == 0, :) = [];
           else
               allEventsTable(eventsIndexTodelete == 1, :) = [];
           end           
       end
    else
        [speedBehave, accelBehave, ~, ~, BehaveDataTreadmil] = treadmilBehave(behaveFileTreadMillPath, globalParameters.behaveFrameRateTM, globalParameters.ImageSamplingRate);
        
        if all(~strcmp(globalParameters.runByEvent, 'non'))
           eventsIndexTodelete = zeros(1, size(allEventsTable, 1));
           
           for i_e = 1:size(allEventsTable, 1)
               
               check_runByEvents = zeros(size(globalParameters.runByEvent));
               for ind = 1:length(globalParameters.runByEvent)
                   if isempty(find(allEventsTable.start(i_e) == BehaveDataTreadmil.(globalParameters.runByEvent{ind}), 1))
                       check_runByEvents(ind) = 1;
                   end
               end
               
               if all(check_runByEvents == 1)
                   eventsIndexTodelete(i_e) = 1;               
               end
               
           end
           
           allEventsTable(eventsIndexTodelete == 1, :) = [];
        else
            
            save([outputpath, 'BehaveTreadMilResults'], 'BehaveDataTreadmil');
        end
        
    end
    
    
    snapnow;
    close all;

    
    if isempty(allEventsTable)
        snapnow;
        close all;
        fclose('all');
        return;
    end
    
    roiActivity_comb = double(roiActivity_comb);
    
%     for HS Activity! with cluster
   
    snapnow;
    close all;
    
    
    for i_cluster = -1:globalParameters.clusterCount
        if i_cluster == 0
            roiActivityPeakSize = 'All';
        elseif i_cluster == -1
            roiActivityPeakSize = 'All_ExcludeBigEvents';
            
        else
            roiActivityPeakSize = ['cluster', num2str(i_cluster)];
            
        end
        
        
        sprintf('Ca Events, %s', roiActivityPeakSize)
            
           
    %     Calc Distance Matrix for ROI in Activity
       switch(globalParameters.roiActivityDistanceFunction)
           case 'WindowEventFULLPearson'
               [roiActivityDistanceMatrixByH, roiActivityDistanceMatrixByP] = calcROIDistanceInActivity_WindowEventPearson_V3(roiActivity_comb, roiActivityNames, selectedROI, allEventsTable, 'FULL', i_cluster, globalParameters.clusterCount);
               
           case 'WindoEventToPeakPearson'
               [roiActivityDistanceMatrixByH, roiActivityDistanceMatrixByP] = calcROIDistanceInActivity_WindowEventPearson_V3(roiActivity_comb, roiActivityNames, selectedROI, allEventsTable, 'ToPeak', i_cluster, globalParameters.clusterCount);

           case 'PeaksPearson'
               [roiActivityDistanceMatrixByH, roiActivityDistanceMatrixByP] = calcROIDistanceInActivity_WindowEventPearson_V3(roiActivity_comb, roiActivityNames, selectedROI, allEventsTable, 'Peaks', i_cluster, globalParameters.clusterCount);
           case 'WindoEventToPeakCov'
               [roiActivityDistanceMatrixByH, roiActivityDistanceMatrixByP] = calcROIDistanceInActivity_WindowEventCov_V3(roiActivity_comb, roiActivityNames, selectedROI, allEventsTable, 'ToPeak', i_cluster, globalParameters.clusterCount);
            
           case 'WindoEventToPeakSperman'
               [roiActivityDistanceMatrixByH, roiActivityDistanceMatrixByP] = calcROIDistanceInActivity_WindowEventSperman_V3(roiActivity_comb, roiActivityNames, selectedROI, allEventsTable, 'ToPeak', i_cluster, globalParameters.clusterCount);
       end
       
       
       sprintf('centrality Analysis By H');
       plotCentralityForGraph(gRoi, roiActivityDistanceMatrixByH, selectedROI, fullfile(centralityFolder, roiActivityPeakSize, 'ByH'), true, selectedROISplitDepth1);

       sprintf('centrality Analysis By P');
       plotCentralityForGraph(gRoi, roiActivityDistanceMatrixByP, selectedROI, fullfile(centralityFolder, roiActivityPeakSize, 'ByP'), true, selectedROISplitDepth1);

       snapnow;
       fclose('all');
       close all;
    end  
    
    snapnow;
    
    close all;
    fclose all;
    
end
