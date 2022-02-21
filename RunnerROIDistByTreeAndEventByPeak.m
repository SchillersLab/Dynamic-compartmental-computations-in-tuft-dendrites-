function RunnerROIDistByTreeAndEventByPeak(neuronTreeFile, activityfile, outputpath)
    [gRoi, ~, selectedROITable] = loadSwcFile(neuronTreeFile, outputpath);
    
    selectedROI = selectedROITable.Name;
    mainTreeBranchROI = selectedROITable.MainBranchID;
    
    
    [roiActivity, roiActivityNames] = loadActivityFile(activityfile, selectedROI);
    
    [distSmallestROiMatrix, clusterSortROI] = calcSmallestDistBetweenRoiInNeuronTree(gRoi, selectedROI, outputpath);

    RunnerLowEventActivity(roiActivity, roiActivityNames, clusterSortROI, distSmallestROiMatrix, selectedROI, mainTreeBranchROI, outputpath, gRoi);
    RunnerHighEventActivity(roiActivity, roiActivityNames, clusterSortROI, distSmallestROiMatrix, selectedROI, mainTreeBranchROI, outputpath, gRoi);    
    RunnerMediumEventActivity(roiActivity, roiActivityNames, clusterSortROI, distSmallestROiMatrix, selectedROI, mainTreeBranchROI, outputpath, gRoi);
    RunnerAllEventActivity(roiActivity, roiActivityNames, clusterSortROI, distSmallestROiMatrix, selectedROI, mainTreeBranchROI, outputpath, gRoi);

end


function RunnerAllEventActivity(roiActivity, roiActivityNames, clusterSortROI, distSmallestROiMatrix, selectedROI, mainTreeBranchROI, outputpath, gRoi)
    
    [eventsWindowsActivity_all, eventsWindowsActivity_toPeak, eventsWindowsActivity_fromPeak, pks_loc] = findEventOfActivity(roiActivity, outputpath);

     plotFunctionsForResults(gRoi, 'PeakEventAll', 'TreeDist', eventsWindowsActivity_toPeak, distSmallestROiMatrix, outputpath, roiActivity, roiActivityNames, selectedROI, pks_loc, clusterSortROI, mainTreeBranchROI);
  end

function RunnerLowEventActivity(roiActivity, roiActivityNames, clusterSortROI, distSmallestROiMatrix, selectedROI, mainTreeBranchROI, outputpath, gRoi)
    [eventsWindowsActivity_all, eventsWindowsActivity_toPeak, eventsWindowsActivity_fromPeak, pks_loc] = findEventOfActivityLow(roiActivity, outputpath);

     plotFunctionsForResults(gRoi, 'PeakEventLow', 'TreeDist', eventsWindowsActivity_toPeak, distSmallestROiMatrix, outputpath, roiActivity, roiActivityNames, selectedROI, pks_loc, clusterSortROI, mainTreeBranchROI);
  end

function RunnerHighEventActivity(roiActivity, roiActivityNames, clusterSortROI, distSmallestROiMatrix, selectedROI, mainTreeBranchROI, outputpath, gRoi)
    
  [eventsWindowsActivity_all, eventsWindowsActivity_toPeak, eventsWindowsActivity_fromPeak, pks_loc] = findEventOfActivityHigh(roiActivity, outputpath);
  plotFunctionsForResults(gRoi, 'PeakEventHigh', 'TreeDist', eventsWindowsActivity_toPeak, distSmallestROiMatrix, outputpath, roiActivity, roiActivityNames, selectedROI, pks_loc, clusterSortROI, mainTreeBranchROI);
  end

function RunnerMediumEventActivity(roiActivity, roiActivityNames, clusterSortROI, distSmallestROiMatrix, selectedROI, mainTreeBranchROI, outputpath, gRoi)    
    [eventsWindowsActivity_all, eventsWindowsActivity_toPeak, eventsWindowsActivity_fromPeak, pks_loc] = findEventOfActivityMedium(roiActivity, outputpath);
     plotFunctionsForResults(gRoi, 'PeakEventMedium', 'TreeDist', eventsWindowsActivity_toPeak, distSmallestROiMatrix, outputpath, roiActivity, roiActivityNames, selectedROI, pks_loc, clusterSortROI, mainTreeBranchROI);
end

function plotFunctionsForResults(gRoi, windowType,distType, eventsWindowsActivity, distROI, outputpath, roiActivity, roiActivityNames, selectedROI, pks_loc, clusterSortROI, mainTreeBranchROI)

      eventsList = calcDistBetweenPeakROIEventsActivity(roiActivity, roiActivityNames, selectedROI, eventsWindowsActivity, outputpath, pks_loc);
      plotDistMatrixForNeuronDistAndAllEventsActivityDist(eventsList, distROI, outputpath, distType, windowType, selectedROI, clusterSortROI, mainTreeBranchROI, gRoi);
    
      corrBetweenNeuronTreeDistAndActivityDistForAllEvents(eventsList, distROI, outputpath, distType, windowType);
end

function [distMatrix, clusterSortROI] = calcSmallestDistBetweenRoiInNeuronTree(gRoi, selectedROI, outputpath)    
    distMatrix = zeros(length(selectedROI), length(selectedROI));
    tickLabels = [];
    
    for index = 1:length(selectedROI)
        tickLabels{index} = selectedROI{index};
        for secIndex = 1:length(selectedROI)
            [~, d] = shortestpath(gRoi,selectedROI{index},selectedROI{secIndex});
            distMatrix(index, secIndex) = d;
        end
    end
    
    clusterNumber = sum(gRoi.Nodes.Depth(: ,1) == 2) * 3;
    
    y = squareform(distMatrix);
    l = linkage(y, 'single');
    c = cluster(l,'maxclust',clusterNumber);
          
    figDendrogram = figure;
    dendrogram(l, 'Labels', tickLabels, 'ColorThreshold', 'default');
    mysave(figDendrogram, [outputpath, '\DendrogramROISmallestDist']);
    
    [~, clusterSortROI] = sort(c);
    
    figDist = figure;
    hold on;
    title({'ROI Distance'});
    xticks(1:length(selectedROI));
    yticks(1:length(selectedROI));
    imagesc(distMatrix(clusterSortROI, clusterSortROI));
    colorbar
    xticklabels(tickLabels(clusterSortROI));
    xtickangle(90);
    yticklabels(tickLabels(clusterSortROI));
    
    mysave(figDist, [outputpath, '\DistMatrixROISmallestDist']);
 
end


function corrBetweenNeuronTreeDistAndActivityDistForAllEvents(eventsList, distROiMatrix, outputpath, type, windowType)
    positiveCorrEvents = 0;
    for events_index = 1:length(eventsList)
        corrMatrix = [];
        corrIndex = 1;
        for index = 1: size(eventsList{events_index}.distMatrix, 2)
            for secIndex = (index + 1): size(eventsList{events_index}.distMatrix, 2)
                corrMatrix(corrIndex, 1:2) = [distROiMatrix(index, secIndex), eventsList{events_index}.distMatrix(index, secIndex)];
                corrIndex = corrIndex + 1;
            end
        end
        
        [RHO,PVAL] = corr(corrMatrix,'type','Pearson');
        if(RHO(1,2) > 0 & PVAL(1,2) <= 0.05)
            positiveCorrEvents = positiveCorrEvents + 1;
        end        
    end
    
    positiveCorrEvents = (positiveCorrEvents / length(eventsList));
    save([outputpath, '\PositiveCorrEventsDistVSDendriticDist_type' type, 'window', windowType], 'positiveCorrEvents');
end


function plotDistMatrixForNeuronDistAndAllEventsActivityDist(eventsList, distROiMatrix, outputpath, type, windowType, selectedROI, clusterSortROI, mainTreeBranchROI, gRoi)
    
    averageEventDistMatrix = zeros(size(distROiMatrix, 1), size(distROiMatrix, 2));
    for events_index = 1:length(eventsList)
        averageEventDistMatrix = averageEventDistMatrix + eventsList{events_index}.distMatrix;
    end
    averageEventDistMatrix = averageEventDistMatrix ./ length(eventsList);
    
    figDist = figure;
    hold on;
    title({'ROI Distance'});
    xticks(1:length(selectedROI));
    yticks(1:length(selectedROI));
    imagesc(averageEventDistMatrix(clusterSortROI, clusterSortROI));
    colorbar
    xticklabels(selectedROI(clusterSortROI));
    xtickangle(90);
    yticklabels(selectedROI(clusterSortROI));
    
    mysave(figDist, [outputpath, '\DistMatrixROIEventsAverageAccordingToTreeDistSort']);
    
    fig = figure;
    hold on;
    % Create ylabel
    ylabel({'Calcium Event Distance'});

    % Create xlabel
    xlabel({'Dendritic distance'});

    title({'Distance Matrix For Average Events'});
    leg = zeros(3, 1);
    leg(1) = plot(0,0, 'color', 'red');
    legColor{1} = 'diffrent main branch';
    
    classesM = unique(mainTreeBranchROI);
    classesColor = {'blue', 'green'};
    
    for index = 1:length(classesM)
        leg(index+1) = plot(0,0, 'color', classesColor{index});
        legColor{index + 1} = gRoi.Nodes(classesM(index),:).Name{1};
    end
    
    for index = 1: size(averageEventDistMatrix, 2)
        for secIndex = (index + 1): size(averageEventDistMatrix, 2)
            color = 'red';
            
            for clr = 1:length(classesM)
                if (mainTreeBranchROI(index) == mainTreeBranchROI(secIndex) && mainTreeBranchROI(secIndex) == classesM(clr))
                    color = classesColor{clr};
                end
            end
            
            hold on;
            scatter(distROiMatrix(index, secIndex), averageEventDistMatrix(index, secIndex), 'filled', 'MarkerFaceColor', color);
        end
    end

    legend(leg, legColor);
    
    mysave(fig, [outputpath, '\EventDistVSDendriticDistForAllEvent_type' type, 'window', windowType]);   
end

function eventsList = calcDistBetweenPeakROIEventsActivity(roiActivity, roiActivityNames, selectedROITreeOrder, eventsWindowsActivity, outputpath, pks_loc)
    for event_index = 1:size(eventsWindowsActivity, 1)
        eventsList{event_index}.window = eventsWindowsActivity(event_index, :);
        eventsList{event_index}.name = ['Event' num2str(event_index)];
        eventsList{event_index}.distMatrix = zeros(length(selectedROITreeOrder), length(selectedROITreeOrder));
    
        for index = 1:length(selectedROITreeOrder)
            activitylocation = strcmpi(roiActivityNames, selectedROITreeOrder{index});
            currentROIActivity = roiActivity(pks_loc(event_index), activitylocation);
              
            for secIndex = 1:length(selectedROITreeOrder)
                activitySeclocation = strcmp(roiActivityNames, selectedROITreeOrder{secIndex});
                secROIActivity = roiActivity(pks_loc(event_index), activitySeclocation);
                
                eventsList{event_index}.distMatrix(index, secIndex) = abs(currentROIActivity(1) - secROIActivity(1));
            end
        end     
   end
end

function [eventsWindowsActivity_all, eventsWindowsActivity_toPeak, eventsWindowsActivity_fromPeak, loc] = findEventOfActivity(roiActivity, outputpath)
    activity_average = mean(roiActivity,2);
    baseLine = iqr(activity_average);
    threshold = 3*baseLine;
    [pks,loc,~,~]  = findpeaks(activity_average,'MinPeakHeight',threshold, 'MinPeakWidth', 3, 'MinPeakDistance', 3);
    
    fig = figure;
    hold on;
    plot(1:length(activity_average), activity_average);  
    plot(loc, pks,'o','MarkerSize',12);
    plot(1:length(activity_average), threshold * ones(length(activity_average),1));
    
    mysave(fig, [outputpath, '\SelectedEventsForROIS']);
  
    eventsWindowsActivity_all = zeros(length(loc), 2);    
    eventsWindowsActivity_toPeak = zeros(length(loc), 2);    
    eventsWindowsActivity_fromPeak = zeros(length(loc), 2);    
    
    for index = 1:length(loc)
        eventsWindowsActivity_all(index, 1) = loc(index) - find(activity_average(loc(index):-1:1) <= baseLine,1) + 1;
        
        belowBaseLine = find(activity_average(loc(index):1:length(activity_average)) <= baseLine,1);
        if isempty(belowBaseLine)
            [~, belowBaseLine] = min(activity_average(loc(index):1:length(activity_average)));
        end
        
        eventsWindowsActivity_all(index, 2) = loc(index) + belowBaseLine - 1;
        
        if (index > 1 && (eventsWindowsActivity_all(index, 1) <= loc(index - 1, 1)))
            act = activity_average(loc(index-1):1:loc(index));
            TF = loc(index-1) + find(act == min(act)) - 1;
            eventsWindowsActivity_all(index, 1) = TF;
            eventsWindowsActivity_all(index - 1, 2) = TF;  
            eventsWindowsActivity_fromPeak(index - 1, 2) = TF;
        end
        
        eventsWindowsActivity_toPeak(index, 1) = eventsWindowsActivity_all(index, 1);
        eventsWindowsActivity_toPeak(index, 2) = loc(index);
        
        eventsWindowsActivity_fromPeak(index, 1) = loc(index);
        eventsWindowsActivity_fromPeak(index, 2) = eventsWindowsActivity_all(index, 2); 
    end
end

function [eventsWindowsActivity_all, eventsWindowsActivity_toPeak, eventsWindowsActivity_fromPeak, loc] = findEventOfActivityLow(roiActivity, outputpath)
    activity_average = mean(roiActivity,2);
    baseLine = iqr(activity_average);
    threshold = 3*baseLine;
    [pks,loc,~,~]  = findpeaks(activity_average,'MinPeakHeight',threshold, 'MinPeakWidth', 3, 'MinPeakDistance', 3);
    
    loc = loc(pks < 2);
    pks = pks(pks < 2);
    
    fig = figure;
    hold on;
    plot(1:length(activity_average), activity_average);  
    plot(loc, pks,'o','MarkerSize',12);
    plot(1:length(activity_average), threshold * ones(length(activity_average),1));
    
    mysave(fig, [outputpath, '\SelectedEventsForROISLow']);
  
    eventsWindowsActivity_all = zeros(length(loc), 2);    
    eventsWindowsActivity_toPeak = zeros(length(loc), 2);    
    eventsWindowsActivity_fromPeak = zeros(length(loc), 2);    
    
    for index = 1:length(loc)
        eventsWindowsActivity_all(index, 1) = loc(index) - find(activity_average(loc(index):-1:1) <= baseLine,1) + 1;
        
        belowBaseLine = find(activity_average(loc(index):1:length(activity_average)) <= baseLine,1);
        if isempty(belowBaseLine)
            [~, belowBaseLine] = min(activity_average(loc(index):1:length(activity_average)));
        end
        
        eventsWindowsActivity_all(index, 2) = loc(index) + belowBaseLine - 1;
        
        if (index > 1 && (eventsWindowsActivity_all(index, 1) <= loc(index - 1, 1)))
            act = activity_average(loc(index-1):1:loc(index));
            TF = loc(index-1) + find(act == min(act)) - 1;
            eventsWindowsActivity_all(index, 1) = TF;
            eventsWindowsActivity_all(index - 1, 2) = TF;  
            eventsWindowsActivity_fromPeak(index - 1, 2) = TF;
        end
        
        eventsWindowsActivity_toPeak(index, 1) = eventsWindowsActivity_all(index, 1);
        eventsWindowsActivity_toPeak(index, 2) = loc(index);
        
        eventsWindowsActivity_fromPeak(index, 1) = loc(index);
        eventsWindowsActivity_fromPeak(index, 2) = eventsWindowsActivity_all(index, 2); 
    end
end

function [eventsWindowsActivity_all, eventsWindowsActivity_toPeak, eventsWindowsActivity_fromPeak, loc] = findEventOfActivityHigh(roiActivity, outputpath)
    activity_average = mean(roiActivity,2);
    baseLine = iqr(activity_average);
    threshold = 3*baseLine;
    [pks,loc,~,~]  = findpeaks(activity_average,'MinPeakHeight',threshold, 'MinPeakWidth', 3, 'MinPeakDistance', 3);
    
    loc = loc(pks > 5);
    pks = pks(pks > 5);
    
    fig = figure;
    hold on;
    plot(1:length(activity_average), activity_average);  
    plot(loc, pks,'o','MarkerSize',12);
    plot(1:length(activity_average), threshold * ones(length(activity_average),1));
    
    mysave(fig, [outputpath, '\SelectedEventsForROISHigh']);
  
    eventsWindowsActivity_all = zeros(length(loc), 2);    
    eventsWindowsActivity_toPeak = zeros(length(loc), 2);    
    eventsWindowsActivity_fromPeak = zeros(length(loc), 2);    
    
    for index = 1:length(loc)
        eventsWindowsActivity_all(index, 1) = loc(index) - find(activity_average(loc(index):-1:1) <= baseLine,1) + 1;
        
        belowBaseLine = find(activity_average(loc(index):1:length(activity_average)) <= baseLine,1);
        if isempty(belowBaseLine)
            [~, belowBaseLine] = min(activity_average(loc(index):1:length(activity_average)));
        end
        
        eventsWindowsActivity_all(index, 2) = loc(index) + belowBaseLine - 1;
        
        if (index > 1 && (eventsWindowsActivity_all(index, 1) <= loc(index - 1, 1)))
            act = activity_average(loc(index-1):1:loc(index));
            TF = loc(index-1) + find(act == min(act)) - 1;
            eventsWindowsActivity_all(index, 1) = TF;
            eventsWindowsActivity_all(index - 1, 2) = TF;  
            eventsWindowsActivity_fromPeak(index - 1, 2) = TF;
        end
        
        eventsWindowsActivity_toPeak(index, 1) = eventsWindowsActivity_all(index, 1);
        eventsWindowsActivity_toPeak(index, 2) = loc(index);
        
        eventsWindowsActivity_fromPeak(index, 1) = loc(index);
        eventsWindowsActivity_fromPeak(index, 2) = eventsWindowsActivity_all(index, 2); 
    end
end

function [eventsWindowsActivity_all, eventsWindowsActivity_toPeak, eventsWindowsActivity_fromPeak, loc] = findEventOfActivityMedium(roiActivity, outputpath)
    activity_average = mean(roiActivity,2);
    baseLine = iqr(activity_average);
    threshold = 3*baseLine;
    [pks,loc,~,~]  = findpeaks(activity_average,'MinPeakHeight',threshold, 'MinPeakWidth', 3, 'MinPeakDistance', 3);
    
    loc = loc(pks >= 2 & pks <= 5);
    pks = pks(pks >= 2 & pks <= 5);
    
    fig = figure;
    hold on;
    plot(1:length(activity_average), activity_average);  
    plot(loc, pks,'o','MarkerSize',12);
    plot(1:length(activity_average), threshold * ones(length(activity_average),1));
    
    mysave(fig, [outputpath, '\SelectedEventsForROISMedium']);
  
    eventsWindowsActivity_all = zeros(length(loc), 2);    
    eventsWindowsActivity_toPeak = zeros(length(loc), 2);    
    eventsWindowsActivity_fromPeak = zeros(length(loc), 2);    
    
    for index = 1:length(loc)
        eventsWindowsActivity_all(index, 1) = loc(index) - find(activity_average(loc(index):-1:1) <= baseLine,1) + 1;
        
        belowBaseLine = find(activity_average(loc(index):1:length(activity_average)) <= baseLine,1);
        if isempty(belowBaseLine)
            [~, belowBaseLine] = min(activity_average(loc(index):1:length(activity_average)));
        end
        
        eventsWindowsActivity_all(index, 2) = loc(index) + belowBaseLine - 1;
        
        if (index > 1 && (eventsWindowsActivity_all(index, 1) <= loc(index - 1, 1)))
            act = activity_average(loc(index-1):1:loc(index));
            TF = loc(index-1) + find(act == min(act)) - 1;
            eventsWindowsActivity_all(index, 1) = TF;
            eventsWindowsActivity_all(index - 1, 2) = TF;  
            eventsWindowsActivity_fromPeak(index - 1, 2) = TF;
        end
        
        eventsWindowsActivity_toPeak(index, 1) = eventsWindowsActivity_all(index, 1);
        eventsWindowsActivity_toPeak(index, 2) = loc(index);
        
        eventsWindowsActivity_fromPeak(index, 1) = loc(index);
        eventsWindowsActivity_fromPeak(index, 2) = eventsWindowsActivity_all(index, 2); 
    end
end