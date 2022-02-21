function [roiActivityDistanceMatrix] = calcROIDistanceInActivity_Euclidean(roiActivity, roiActivityNames, selectedROI, windowFULL)
    for event_index = 1:size(windowFULL, 1)
        eventsList{event_index}.window = windowFULL(event_index, :);
        eventsList{event_index}.name = ['Event' num2str(event_index)];
        eventsList{event_index}.distMatrix = zeros(length(selectedROI), length(selectedROI));
    
        windowIndexArray = eventsList{event_index}.window(1):eventsList{event_index}.window(2);    
        for index = 1:length(selectedROI)
            activitylocation = strcmpi(roiActivityNames, selectedROI{index});
            currentROIActivity = roiActivity(windowIndexArray, activitylocation);
                        
            for secIndex = 1:length(selectedROI)
                activitySeclocation = strcmp(roiActivityNames, selectedROI{secIndex});
                secROIActivity = roiActivity(windowIndexArray, activitySeclocation);
                
                eventsList{event_index}.distMatrix(index, secIndex) = norm(currentROIActivity - secROIActivity); 
            end
        end
    end
    
    roiActivityDistanceMatrix = zeros(size(selectedROI, 1), size(selectedROI, 1));
    for events_index = 1:length(eventsList)
        roiActivityDistanceMatrix = roiActivityDistanceMatrix + eventsList{events_index}.distMatrix;
    end
    roiActivityDistanceMatrix = roiActivityDistanceMatrix ./ length(eventsList);
   
end