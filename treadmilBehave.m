function [speedResultsBehave, accelResultsBehave, speedResultsActivity, accelResultsActivity, BehaveDataTreadmil] = treadmilBehave(treadmilTxtFile, behaveFrameRate, activityFrameR)
    treadmilData = readtable(treadmilTxtFile);
    
    treadmilDataNormal = createTreadMilData(behaveFrameRate, treadmilData);
    treadmilDataEvents = createTreadMilData(behaveFrameRate*2, treadmilData);
    
    treadmilDataNormal.twoP = seconds(treadmilDataNormal.twoP);
    treadmilASTimeTable = table2timetable(treadmilDataNormal, 'RowTimes', 'twoP');
    treadmilASTimeTable = retime(treadmilASTimeTable, 'secondly', 'mean');
    
%     treadmilDataEvents.twoP = seconds(treadmilDataEvents.twoP);
%     treadmilASTimeTableE = table2timetable(treadmilDataEvents, 'RowTimes', 'twoP');
%     treadmilASTimeTableE = retime(treadmilASTimeTableE, 'secondly', 'mean');
%     
    treadmilASTimeTableE = treadmilASTimeTable;

    speedResultsBehave = treadmilASTimeTable.speed;
    accelResultsBehave = treadmilASTimeTable.accel;
    
    speedResultsActivity = treadmilDataNormal.speed;
    accelResultsActivity = treadmilDataNormal.accel;
    
    speedResultsActivityEvents = treadmilASTimeTableE.speed;
    accelResultsActivityEvents = treadmilASTimeTableE.accel;    
    
    std_forbaseline_velocity = std(speedResultsActivityEvents(~isnan(speedResultsActivityEvents)& speedResultsActivityEvents < 0.1));
    speedResultsActivityEvents(speedResultsActivityEvents < std_forbaseline_velocity & speedResultsActivityEvents > -1*std_forbaseline_velocity) = 0;
    
    std_forbaseline_Acc = std(accelResultsActivityEvents(~isnan(accelResultsActivityEvents)& accelResultsActivityEvents < 0.1));
    accelResultsActivityEvents(accelResultsActivityEvents < std_forbaseline_Acc & accelResultsActivityEvents > -1*std_forbaseline_Acc) = 0;
%     
    slop_build_velocity = zeros(size(speedResultsActivityEvents));
    slop_build_velocity(speedResultsActivityEvents ~= 0) = 1;
    
    slop_build_velocity = slopFix(slop_build_velocity, speedResultsActivityEvents, activityFrameR);
    
    slop_build_acc = zeros(size(accelResultsActivityEvents));
    slop_build_acc(accelResultsActivityEvents > 0) = 1;
    slop_build_acc(accelResultsActivityEvents < 0) = -1;
    
%     slop_build_acc = slopFix(slop_build_acc, accelResultsActivityEvents, activityFrameR);
    
    BehaveDataTreadmil.rest = find(slop_build_velocity == 0 & slop_build_acc == 0);
    BehaveDataTreadmil.walkconstant = find(slop_build_velocity ~= 0 & slop_build_acc == 0);
    BehaveDataTreadmil.walkacceleration = find(slop_build_velocity ~= 0 & slop_build_acc ~= 0);
    BehaveDataTreadmil.walkPosacceleration = find(slop_build_velocity ~= 0 & slop_build_acc > 0);
    BehaveDataTreadmil.walkNegacceleration = find(slop_build_velocity ~= 0 & slop_build_acc < 0);
    
    res = find(slop_build_velocity(1:end-1) ~= slop_build_velocity(2:end));
    BehaveDataTreadmil.stop = [];
    BehaveDataTreadmil.start = [];
    for ind = 1:length(res)
        if slop_build_velocity(res(ind)) == 1
           stop_vector = max(1, res(ind) - 4) : min(length(slop_build_velocity), res(ind) + 5);
           BehaveDataTreadmil.stop = [BehaveDataTreadmil.stop, stop_vector];
        else
            start_vector = max(1, res(ind) - 4) : min(length(slop_build_velocity), res(ind) + 5);
            BehaveDataTreadmil.start = [BehaveDataTreadmil.start, start_vector];
        end
    end    
end

function treadmilData = createTreadMilData(behaveFrameRate, treadmilData)
    
    b_floor = floor(behaveFrameRate / 2);
    b_celing = ceil(behaveFrameRate / 2);
    
    deltaTreadmilSpeed = (treadmilData.treadmill((1 + behaveFrameRate):size(treadmilData,1)) - treadmilData.treadmill(1:(size(treadmilData,1)-behaveFrameRate)));
    speed_vector =  (deltaTreadmilSpeed / 1024) * pi * 10;
    treadmilData.speed = [NaN(b_floor, 1); speed_vector; NaN(b_celing, 1)];
    
    deltaTreadmilAccel = (treadmilData.speed((1 + behaveFrameRate):size(treadmilData,1)) - treadmilData.speed(1:(size(treadmilData,1)-behaveFrameRate)));
    
    treadmilData.accel = [NaN(b_floor,1); deltaTreadmilAccel; NaN(b_celing,1)];
    
    [maxTwoP_v, maxTwoP_i] = max(treadmilData.twoP);
    
    zeroAfterMax = find(treadmilData.twoP((maxTwoP_i + 1) : end) < maxTwoP_v) + maxTwoP_i;
    treadmilData(zeroAfterMax, :) = [];
end

function slop_build = slopFix(slop_build, speedResultsActivity, activityFrameR)
    res = [0; find(slop_build(1:end-1) ~= slop_build(2:end)); length(speedResultsActivity)];
    
    for k = 1:length(res)-1
       if slop_build(res(k) + 1) == 0 & (res(k + 1) - res(k) <= 3)
            if (k-1>0)
                slop_value = slop_build(res(k-1) + 1) ;    
            else
                slop_value = slop_build(res(k+1) + 1) ;
            end
            
            slop_build((res(k) + 1):res(k+1)) = slop_value;
        end
    end
   
    for k = 1:length(res)-1
       if (slop_build(res(k) + 1) == 1 | slop_build(res(k) + 1) == -1) & (res(k + 1) - res(k) <= 3)
            slop_value = 0;
            if (k-1>0)
                slop_value = slop_build(res(k-1) + 1) ;    
            end
            
            slop_build((res(k) + 1):res(k+1)) = slop_value;
        end
    end
    
end