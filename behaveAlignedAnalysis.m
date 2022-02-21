function behaveAlignedAnalysis(globalParameters)
    eventsDetectionFolder = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, ...
        globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
        globalParameters.RunnerDate,globalParameters.RunnerNumber, 'EventsDetection');
    load([eventsDetectionFolder, '\roiActivity_comb.mat'], 'allEventsTable', 'roiActivity_comb'); 
    neuronActivityPathTPA = fullfile(globalParameters.TPAFolder, globalParameters.AnimalName, globalParameters.DateAnimal);
    
    outputpath = fullfile('\\jackie-analysis\e\Shay\StatisticSummary\GLMBehaveAligned', [globalParameters.AnimalName globalParameters.DateAnimal globalParameters.neuronNumberName]);
    mkdir(outputpath);
    
    cluster1SumOriginal = sum(allEventsTable.clusterByH == 1);
    cluster2SumOriginal = sum(allEventsTable.clusterByH == 2);
    cluster3SumOriginal = sum(allEventsTable.clusterByH == 3);
    cluster4SumOriginal = sum(allEventsTable.clusterByH == 4);
    
    tr_frame_count = globalParameters.ImageSamplingRate * globalParameters.time_sec_of_trial;
    
    [BehaveDataAll, NAMES, trials_label] = loadBDAFile(neuronActivityPathTPA, globalParameters.BehavioralSamplingRate, globalParameters.ImageSamplingRate, tr_frame_count, globalParameters.behavioralDelay, globalParameters.toneTime);
    
    for i_e = 1:size(allEventsTable, 1)
        tr_index = floor(allEventsTable.start(i_e) ./ tr_frame_count) + 1;
        allEventsTable.tr_index(i_e) = tr_index;
    end
    
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
      plotEventCaForBehaveDataHandReach(BehaveDataAll.(NAMES{currentEventLoc(in)}).startTiming, tr_frame_count, allEventsTable, globalParameters.clusterCount, outputpath, NAMES{currentEventLoc(in)}, trials_label, globalParameters.split_trialsLabel)           
   end
   
   allEventsTable.behave(:) = cell(1, size(allEventsTable, 1));
          
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

       allEventsTable.behave(i_e) = {runByEventTemp(checkEvent ~= 1)};


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
   
   cluster1Precentage = sum(allEventsTable.clusterByH == 1) ./ cluster1SumOriginal;   
   cluster2Precentage = sum(allEventsTable.clusterByH == 2) ./ cluster2SumOriginal;
   cluster3Precentage = sum(allEventsTable.clusterByH == 3) ./ cluster3SumOriginal;
   cluster4Precentage = sum(allEventsTable.clusterByH == 4) ./ cluster4SumOriginal;
   
   save(fullfile('\\jackie-analysis\e\Shay\StatisticSummary\GLMBehaveAligned', [globalParameters.AnimalName '_' globalParameters.DateAnimal '_' globalParameters.neuronNumberName '.mat']), 'cluster1Precentage', 'cluster2Precentage', 'cluster3Precentage', 'cluster4Precentage');
end