function [roiActivityDistanceMatrix] = calcROIDistanceInActivity_WindowEventPearson(roiActivity, roiActivityNames, selectedROI, all_event_struct, typeCom, clusterNum)
    for index = 1:length(selectedROI)
        activitylocation = strcmpi(roiActivityNames, selectedROI{index});
                
        locationToCompare = getLocTocompare(activitylocation, typeCom, clusterNum, all_event_struct);
        
        for secIndex = 1:length(selectedROI)
            activitySeclocation = strcmp(roiActivityNames, selectedROI{secIndex});
            locationToCompareSec = getLocTocompare(activitySeclocation, typeCom, clusterNum, all_event_struct);
        
            totalActivityLocation = [locationToCompare, locationToCompareSec];
            totalActivityLocation = sort(totalActivityLocation);
            totalActivityLocation = unique(totalActivityLocation);
            
            currentROIActivity = roiActivity(totalActivityLocation, activitylocation);
            secROIActivity = roiActivity(totalActivityLocation, activitySeclocation);
            
            if isempty(totalActivityLocation)
                roiActivityDistanceMatrix(index, secIndex) = 0;
                continue;
            end
            
            corrEventsPeaksROI = corr([currentROIActivity, secROIActivity], 'type', 'Pearson');
%             if index == secIndex
%                 roiActivityDistanceMatrix(index, secIndex) = 0;
%             else
%                 roiActivityDistanceMatrix(index, secIndex) = 1 - corrEventsPeaksROI(1, 2);
%             end

            roiActivityDistanceMatrix(index, secIndex) = corrEventsPeaksROI(1, 2);
        end
    end
end

function locationToCompare = getLocTocompare(activityLoc, typeCom, clusterNum, all_event_struct)
    switch typeCom
        case 'FULL'
            startEventList = all_event_struct.start{activityLoc};            
            endEventList = all_event_struct.end{activityLoc};
        case 'ToPeak'
            startEventList = all_event_struct.start{activityLoc};            
            endEventList = all_event_struct.pks{activityLoc};
        case 'Peaks'
            startEventList = all_event_struct.pks{activityLoc};            
            endEventList = all_event_struct.pks{activityLoc};
    end
    
    indexByCluster = all_event_struct.cluster{activityLoc};
    
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