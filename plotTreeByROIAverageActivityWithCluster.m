function plotTreeByROIAverageActivityWithCluster(gRoi, outputpath, roiActivity, roiActivityNames, selectedROI, all_event_table, typeCom, clusterNum, roiSorted)
    averageActivityForROIByH = zeros(length(selectedROI), 1);
    averageActivityForROIByP = zeros(length(selectedROI), 1);
    
    for index = 1:length(selectedROI)
        activitylocation = strcmpi(roiActivityNames, selectedROI{index});
        
        locationToCompareByH = getLocTocompare(typeCom, clusterNum, all_event_table,(all_event_table.clusterByH), activitylocation);
        locationToCompareByPrecantage = getLocTocompare(typeCom, clusterNum, all_event_table,(all_event_table.clusterByRoiPrecantage), activitylocation);
        
        currentROIActivityByH = roiActivity(locationToCompareByH, activitylocation);
        currentROIActivityByP = roiActivity(locationToCompareByPrecantage, activitylocation);
        
        if isempty(currentROIActivityByP)
            averageActivityForROIByP(index) = 0;
        else
            averageActivityForROIByP(index) = mean(currentROIActivityByP);
        end
        
        if isempty(currentROIActivityByH)
            averageActivityForROIByH(index) = 0;
        else
            averageActivityForROIByH(index) = mean(currentROIActivityByH);
        end
    end
    
    
    nodesColor = zeros(length(gRoi.Nodes.Name),3);
    
    rgbColors = vals2colormap(averageActivityForROIByH, 'jet', [min(averageActivityForROIByH), max(averageActivityForROIByH)]);
    
    for index = 1:length(selectedROI)
        locRoi = find(strcmp(gRoi.Nodes.Name, selectedROI(index)));
        nodesColor(locRoi, :) = rgbColors(index,:);
    end
     
    fileName2 = [outputpath, '\GraphWithROI_ColoredActivity' 'ByH'];
    
    plotGraphWithROI(gRoi, fileName2, nodesColor, {'Roi ByActivity cluster By H '});
    
    nodesColor = zeros(length(gRoi.Nodes.Name),3);
    
    rgbColors = vals2colormap(averageActivityForROIByP, 'jet', [min(averageActivityForROIByP), max(averageActivityForROIByP)]);
    
    for index = 1:length(selectedROI)
        locRoi = find(strcmp(gRoi.Nodes.Name, selectedROI(index)));
        nodesColor(locRoi, :) = rgbColors(index,:);
    end
     
    
    fileName2 = [outputpath, '\GraphWithROI_ColoredActivity' 'ByP'];
    plotGraphWithROI(gRoi, fileName2, nodesColor, {'Roi ByActivity cluster By P '});
   
    fig = figure;
    plot(averageActivityForROIByH(roiSorted), '-*k');
    xticks(1:length(selectedROI));
    xtickangle(90);
    xticklabels(selectedROI(roiSorted));
    ylabel('mean df/f for cluster events by Hight');
    mysave(fig, [outputpath, '\MeanActivityROI_' 'ByH']);
    
    fig = figure;
    plot(averageActivityForROIByP(roiSorted), '-*k');
    xticks(1:length(selectedROI));
    xtickangle(90);
    xticklabels(selectedROI(roiSorted));
    ylabel('mean df/f for cluster events by Precentage');
    mysave(fig, [outputpath, '\MeanActivityROI_' 'ByP']);
end

function locationToCompare = getLocTocompare(typeCom, clusterNum, all_event_table, indexByCluster, activitylocation)
    eventsIncluded = zeros(1, size(all_event_table, 1));
    for indexEvent = 1:size(all_event_table, 1)
        currentEventROIIndex = all_event_table.roisEvent{indexEvent};
        if sum(activitylocation & currentEventROIIndex') > 0
            eventsIncluded(indexEvent) = 1;
        end
    end

    eventsLocation = eventsIncluded == 1;
    indexByCluster = indexByCluster(eventsLocation);
    
    switch typeCom
        case 'FULL'
            startEventList = all_event_table.start(eventsLocation);            
            endEventList = all_event_table.event_end(eventsLocation);
        case 'ToPeak'
            startEventList = all_event_table.start(eventsLocation);            
            endEventList = all_event_table.pks(eventsLocation);
        case 'Peaks'
            startEventList = all_event_table.pks(eventsLocation);            
            endEventList = all_event_table.pks(eventsLocation);
    end
    
    if clusterNum == 0
        startEventListBycluster = startEventList;
        endEventListBycluster = endEventList;
    elseif clusterNum == -1
        startEventListBycluster = startEventList(indexByCluster ~= 1);
        endEventListBycluster = endEventList(indexByCluster ~= 1);
    else
        startEventListBycluster = startEventList(indexByCluster == clusterNum);
        endEventListBycluster = endEventList(indexByCluster == clusterNum);
    end
    
    locationToCompare = [];
    
    for i = 1:length(startEventListBycluster)
        locationToCompare = [locationToCompare, startEventListBycluster(i):endEventListBycluster(i)];
    end
end