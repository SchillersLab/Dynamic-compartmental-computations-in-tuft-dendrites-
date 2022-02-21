function [startTiming, endTiming] = findFirstEventsAfter(BehaveData, NAMES, eventName, afterName)   
    locationN = find(strcmp(NAMES, eventName));
    
    startTiming = zeros(1, length(BehaveData.(NAMES{locationN}).startTiming));
    endTiming = zeros(1, length(BehaveData.(NAMES{locationN}).startTiming));
    
    for k = 1:length(BehaveData.(NAMES{locationN}).startTiming)
        
        if BehaveData.(NAMES{locationN}).startTiming(k) ~= 0
            eN = findFirstEventsAfterFun(BehaveData, NAMES, BehaveData.(NAMES{locationN}).startTiming(k), k);
            if contains(eN, afterName)
                startTiming(k) = BehaveData.(NAMES{locationN}).startTiming(k);
                endTiming(k) = BehaveData.(NAMES{locationN}).endTiming(k);
            end
        end
    end

end

function eventName = findFirstEventsAfterFun(BehaveData, NAMES, currentEventStartTime, currentEventTrial)
    lastLocation = contains(NAMES, 'last');
    NAMES(lastLocation) = [];
    ftLocation = contains(NAMES, 'firstTone');
    NAMES(ftLocation) = [];
    eventName = 'non';
    lagEventMin = inf;
    for i  = 1:length(NAMES)
        if BehaveData.(NAMES{i}).startTiming(currentEventTrial) ~= 0
            curTimeLag = BehaveData.(NAMES{i}).startTiming(currentEventTrial) - currentEventStartTime;

            if curTimeLag > 0 && curTimeLag < lagEventMin
                lagEventMin = curTimeLag;
                eventName = NAMES{i};
            end
        end
    end
end