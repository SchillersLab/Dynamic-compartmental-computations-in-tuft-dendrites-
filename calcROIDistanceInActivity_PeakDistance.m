function [roiActivityDistanceMatrix] = calcROIDistanceInActivity_PeakDistance(roiActivity, roiActivityNames, selectedROI, locationPeaks)
   for event_index = 1:size(locationPeaks, 1)
        eventsList{event_index}.name = ['Event' num2str(event_index)];
        eventsList{event_index}.distMatrix = zeros(length(selectedROI), length(selectedROI));
    
        for index = 1:length(selectedROI)
            activitylocation = strcmpi(roiActivityNames, selectedROI{index});
            currentROIActivity = roiActivity(locationPeaks(event_index), activitylocation);
              
            for secIndex = 1:length(selectedROI)
                activitySeclocation = strcmp(roiActivityNames, selectedROI{secIndex});
                secROIActivity = roiActivity(locationPeaks(event_index), activitySeclocation);
                
                eventsList{event_index}.distMatrix(index, secIndex) = abs(currentROIActivity(1) - secROIActivity(1));
            end
        end     
   end
   
   roiActivityDistanceMatrix = zeros(size(selectedROI, 1), size(selectedROI, 1));
    for events_index = 1:length(eventsList)
        roiActivityDistanceMatrix = roiActivityDistanceMatrix + eventsList{events_index}.distMatrix;
    end
    roiActivityDistanceMatrix = roiActivityDistanceMatrix ./ length(eventsList);
   
end