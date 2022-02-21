function [speedResultsBehave, accelResultsBehave, speedResultsActivity, accelResultsActivity, BehaveDataTreadmil] = treadmilBehaveType2(treadmilTxtFile, behaveFrameRate, activityFrameR, std_treadMilThreshold, winlengthsec, activityTotalSize)
    treadmilData = readtable(treadmilTxtFile);
    
    maxTwoP = max(treadmilData.twoP);
    lastTwoP = find(treadmilData.twoP == maxTwoP, 1, 'last');
    BehaveSamplingRate = (treadmilData.time(lastTwoP) - treadmilData.time(1)) ./ 1000;
    BehaveSamplingRate = ceil(lastTwoP ./ BehaveSamplingRate);

    
    treadmilDataNormal = createTreadMilData(behaveFrameRate, treadmilData);
    
    treadmilDataNormal.twoP = seconds(treadmilDataNormal.twoP);
    treadmilASTimeTable = table2timetable(treadmilDataNormal, 'RowTimes', 'twoP');
    treadmilASTimeTable = retime(treadmilASTimeTable, 'secondly', 'mean');
    
    treadmilASTimeTable(activityTotalSize+1:end,:) = [];
    
    nanValues = find(isnan(treadmilASTimeTable.treadmill));
    if ~isempty(nanValues)
        nanValuesM = nanValues(2:end) - nanValues(1:end-1);
    
        if all((nanValuesM == 1))
            treadmilASTimeTable.treadmill(nanValues) = treadmilASTimeTable.treadmill(nanValues(end) + 1);
        else
            nanValuesM2 = nanValues([1; find(nanValuesM ~= 1) + 1]);
            treadmilASTimeTable.treadmill(nanValuesM2) = treadmilASTimeTable.treadmill(nanValuesM2 - 1);    
        end    
    end
    
    smooth_speed = smoothdata(treadmilASTimeTable.speed,'gaussian', activityFrameR , 'omitnan');
    Y=nan(size(smooth_speed));

    for i = 1:length(smooth_speed)-winlengthsec*activityFrameR
        Y(i+winlengthsec*activityFrameR/2) = std(smooth_speed(i:(i+winlengthsec*activityFrameR-1)), 'omitnan');
    end
    
    Y1=nan(size(Y));
    Y2=nan(size(Y));
    for i = 1:length(Y)
        Y1(i+winlengthsec*activityFrameR/2) = Y(min(i + winlengthsec*activityFrameR, length(Y))) ./ Y(i);
        Y2(i+winlengthsec*activityFrameR/2) = Y(max(i - winlengthsec*activityFrameR, 1)) ./ Y(i);
    end
       
    Z = Y > std_treadMilThreshold;
    
    runOnset = find(Z(2:end)-Z(1:(end-1)) == 1);
    runOffset = find(Z(2:end)-Z(1:(end-1)) == -1);
    
    runOnset = runOnset + winlengthsec*activityFrameR/2;
    runOffset = runOffset + winlengthsec*activityFrameR/2;
% %     
    isTrueRunSet = zeros(1, length(runOnset));
    BehaveDataTreadmil.walk = [];
    
    for i = 1:length(runOnset)
        [~, maxT] = max(Y1(max(runOnset(i) -  winlengthsec*activityFrameR ./ 2, 0):min(length(Y1), runOnset(i) +  winlengthsec*activityFrameR ./ 2)));
        runOnset(i) = maxT + max(runOnset(i) -  winlengthsec*activityFrameR ./ 2, 0) - 1;
        
        if i > length(runOffset)
            runOffset(i) = length(Z);
        end
        
        [~, minT] = min(Y2(max(runOffset(i) -  winlengthsec*activityFrameR./2, 0):min(length(Y2), runOffset(i) +  winlengthsec*activityFrameR./2)));
        runOffset(i) = minT + max(runOffset(i) -  winlengthsec*activityFrameR./2, 0) - 1;
        
        if (runOffset(i) - runOnset(i) >= activityFrameR)
%         if (runOffset(i) - runOnset(i) >= activityFrameR && nanmean(smooth_speed(runOnset(i):runOffset(i))) >= 0.005)
            isTrueRunSet(i) = 1;
            
            BehaveDataTreadmil.walk(end+1:(end+1+runOffset(i) - runOnset(i))) = runOnset(i):runOffset(i);
        end
    end
    
    runOnset(isTrueRunSet == 0) = [];
    runOffset(isTrueRunSet == 0) = [];
    
    accelSmooth = diff(smooth_speed);
    accelSmooth(end+1) = accelSmooth(end);
    speedResultsBehave = smooth_speed;
    accelResultsBehave = accelSmooth;
    
    speedResultsActivity = treadmilDataNormal.speed;
    accelResultsActivity = treadmilDataNormal.accel;
    
    BehaveDataTreadmil.rest = 1:length(smooth_speed);
    BehaveDataTreadmil.rest(BehaveDataTreadmil.walk) = [];
    
    walkVector = zeros(length(smooth_speed), 1);
    walkVector(BehaveDataTreadmil.walk) = 1;
    
    BehaveDataTreadmil.walkconstant = find(accelSmooth == 0 & walkVector == 1);
    BehaveDataTreadmil.walkacceleration = find(accelSmooth ~= 0 & walkVector == 1);
    BehaveDataTreadmil.walkPosacceleration = find(accelSmooth > 0 & walkVector == 1);
    BehaveDataTreadmil.walkNegacceleration = find(accelSmooth < 0 & walkVector == 1);
    BehaveDataTreadmil.stop = runOffset;
    BehaveDataTreadmil.start = runOnset;
    BehaveDataTreadmil.speed = smooth_speed;
    BehaveDataTreadmil.accel = accelSmooth;
    BehaveDataTreadmil.walkBinary = walkVector;
    BehaveDataTreadmil.restBinary = zeros(length(smooth_speed), 1);
    BehaveDataTreadmil.restBinary(BehaveDataTreadmil.rest) = 1;
    BehaveDataTreadmil.posaccelBinary = double(accelSmooth > 0 & walkVector == 1);
    BehaveDataTreadmil.negaccelBinary = double(accelSmooth < 0 & walkVector == 1);
    BehaveDataTreadmil.x_location = treadmilASTimeTable.treadmill;
    BehaveDataTreadmil.BehaveSamplingRate = BehaveSamplingRate;
end

function treadmilData = createTreadMilData(behaveFrameRate, treadmilData)
    
    b_floor = floor(behaveFrameRate / 2);
    b_celing = ceil(behaveFrameRate / 2);
    
    deltaTreadmilSpeed = (treadmilData.treadmill((1 + behaveFrameRate):size(treadmilData,1)) - treadmilData.treadmill(1:(size(treadmilData,1)-behaveFrameRate))) ./ behaveFrameRate;
    speed_vector =  (deltaTreadmilSpeed / 1024) * pi * 10;
    treadmilData.speed = [NaN(b_floor, 1); speed_vector; NaN(b_celing, 1)];
    
    deltaTreadmilAccel = (treadmilData.speed((1 + behaveFrameRate):size(treadmilData,1)) - treadmilData.speed(1:(size(treadmilData,1)-behaveFrameRate))) ./ behaveFrameRate;
    
    treadmilData.accel = [NaN(b_floor,1); deltaTreadmilAccel; NaN(b_celing,1)];
    
    [maxTwoP_v, maxTwoP_i] = max(treadmilData.twoP);
    
    zeroAfterMax = find(treadmilData.twoP((maxTwoP_i + 1) : end) < maxTwoP_v) + maxTwoP_i;
    treadmilData(zeroAfterMax, :) = [];
end
