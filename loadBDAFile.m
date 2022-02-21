function [BehaveData, NAMES, trials_label] = loadBDAFile(activityfileTPAFolder, BehavioralSamplingRate, ImagingSamplingRate, tr_frame_count, BehavioralDelay, toneTime)
    frameRateRatio = BehavioralSamplingRate/ImagingSamplingRate;
    
    BdaList = dir(strcat(activityfileTPAFolder,'\BDA*'));

    fileNumRoi = length(BdaList);
   
    toneTimeFrame = round((toneTime*BehavioralSamplingRate-BehavioralDelay)./frameRateRatio); % transfers to time of the two photon
    
    eventNameList = [];
    allTrialEvents                = cell(fileNumRoi,1);
    for trialInd = 1:fileNumRoi
        usrData                    = load(fullfile(BdaList(trialInd).folder, BdaList(trialInd).name));
        allTrialEvents{trialInd}   = usrData.strEvent;
        for event_i = 1:length(allTrialEvents{trialInd})
            if isempty(eventNameList) || ~any(strcmpi(eventNameList, allTrialEvents{trialInd}{event_i}.Name))
                eventNameList{end+1} = extractEventstr(allTrialEvents{trialInd}{event_i}.Name);
            end
        end
    end
    
    framNum = tr_frame_count;
    for eventName_i = 1:length(eventNameList)
        BehaveData.(eventNameList{eventName_i}).indicator = zeros(fileNumRoi, framNum); 
        BehaveData.(eventNameList{eventName_i}).startTiming = zeros(fileNumRoi, 1);
        BehaveData.(eventNameList{eventName_i}).endTiming = zeros(fileNumRoi, 1);
        
        if ~contains(eventNameList{eventName_i}, 'success') && ~contains(eventNameList{eventName_i}, 'failure')
            eventTName = ['last' eventNameList{eventName_i}(isletter(eventNameList{eventName_i}))];
            BehaveData.(eventTName).indicator = zeros(fileNumRoi, framNum); 
            BehaveData.(eventTName).startTiming = zeros(fileNumRoi, 1);
            BehaveData.(eventTName).endTiming = zeros(fileNumRoi, 1);
            BehaveData.(eventTName).count = zeros(fileNumRoi, 1);

            eventTName = ['firstTone' eventNameList{eventName_i}(isletter(eventNameList{eventName_i}))];
            BehaveData.(eventTName).indicator = zeros(fileNumRoi, framNum); 
            BehaveData.(eventTName).startTiming = zeros(fileNumRoi, 1);
            BehaveData.(eventTName).endTiming = zeros(fileNumRoi, 1);
            BehaveData.(eventTName).count = zeros(fileNumRoi, 1);
            
             eventTName = ['secTone' eventNameList{eventName_i}(isletter(eventNameList{eventName_i}))];
            BehaveData.(eventTName).indicator = zeros(fileNumRoi, framNum); 
            BehaveData.(eventTName).startTiming = zeros(fileNumRoi, 1);
            BehaveData.(eventTName).endTiming = zeros(fileNumRoi, 1);
            BehaveData.(eventTName).count = zeros(fileNumRoi, 1);
        end
            
    end
    
    NAMES = fieldnames(BehaveData);
        
    for trial_i = 1:fileNumRoi
        for m = 1:length(allTrialEvents{trial_i})
            eventname = lower(allTrialEvents{trial_i}{m}.Name);
            eventname = extractEventstr(eventname);
            if length(allTrialEvents{trial_i}{m}.tInd) ==2
                timeInd     = allTrialEvents{trial_i}{m}.tInd;
            else
                timeInd     = allTrialEvents{trial_i}{m}.TimeInd;
            end
            if isempty(timeInd)
                continue;
            end
            %             frameRateRatio=size(allTrialEvents{trial_i}{end}.Data,1)/size(eventDataArray,1);
            %                 frameRateRatio=18
            timeInd     = round((timeInd-BehavioralDelay)./frameRateRatio); % transfers to time of the two photon
            timeInd     = max(1,min(framNum,timeInd));
            % assign to vector
            
            if contains(eventname, 'success') || contains(eventname, 'failure')
                BehaveData.(eventname).indicator(trial_i, :) = 1;
                BehaveData.(eventname).startTiming(trial_i) = timeInd(1);
                BehaveData.(eventname).endTiming(trial_i) = timeInd(2);
            else
                BehaveData.(eventname).indicator(trial_i, timeInd(1):timeInd(2)) = 1;
                BehaveData.(eventname).startTiming(trial_i) = timeInd(1);
                BehaveData.(eventname).endTiming(trial_i) = timeInd(2);
                
                eventTName = ['last' eventname(isletter(eventname))];

                if BehaveData.(eventTName).startTiming(trial_i) < timeInd(1)
                    BehaveData.(eventTName).indicator(trial_i, :) = 0;
                    BehaveData.(eventTName).indicator(trial_i, timeInd(1):timeInd(2)) = 1;
                    BehaveData.(eventTName).startTiming(trial_i) = timeInd(1);
                    BehaveData.(eventTName).endTiming(trial_i) = timeInd(2);
                    BehaveData.(eventTName).count(trial_i) = BehaveData.(eventTName).count(trial_i) + 1;
                end
                
                eventTName = ['firstTone' eventname(isletter(eventname))];

                if BehaveData.(eventTName).startTiming(trial_i) == 0 && timeInd(1) >= toneTimeFrame
                    BehaveData.(eventTName).indicator(trial_i, :) = 0;
                    BehaveData.(eventTName).indicator(trial_i, timeInd(1):timeInd(2)) = 1;
                    BehaveData.(eventTName).startTiming(trial_i) = timeInd(1);
                    BehaveData.(eventTName).endTiming(trial_i) = timeInd(2);
                    BehaveData.(eventTName).count(trial_i) = BehaveData.(eventTName).count(trial_i) + 1;
                end
                
                eventTName = ['secTone' eventname(isletter(eventname))];

                if BehaveData.(eventTName).startTiming(trial_i) == 0 && timeInd(1) >= toneTimeFrame &&...
                        BehaveData.(['firstTone' eventname(isletter(eventname))]).startTiming(trial_i) < timeInd(1)
                    BehaveData.(eventTName).indicator(trial_i, :) = 0;
                    BehaveData.(eventTName).indicator(trial_i, timeInd(1):timeInd(2)) = 1;
                    BehaveData.(eventTName).startTiming(trial_i) = timeInd(1);
                    BehaveData.(eventTName).endTiming(trial_i) = timeInd(2);
                    BehaveData.(eventTName).count(trial_i) = BehaveData.(eventTName).count(trial_i) + 1;
                end
            end
                       
        end        
    end
    
    trials_label = zeros(1, fileNumRoi);
    trials_label(sum(BehaveData.success.indicator, 2) > 0) = 1;
    trials_label(sum(BehaveData.failure.indicator, 2) > 0) = 2;
end