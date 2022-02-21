function saveDataFotGLMTreadMil_type3(doBehave, neuronNumberName, splinesLocation ,BehaveDataTreadmil, roiActivity, roiActivityNames, ImageSamplingRate, tpaFolder, aftertonetime, treadMilTxt, videoPath)
    treadmilData = readtable(treadMilTxt);   
    
    [filepath,name,ext] = fileparts(tpaFolder);
    
    outputpath = fullfile(filepath, [name, ext, '_glm']);
    
    videoTimeInsec = ((treadmilData.time(end) - treadmilData.time(1)) ./  1000);
    videoSamplingRate = ceil(size(treadmilData, 1) ./ videoTimeInsec);
    
    [maxTwoP_v, maxTwoP_i] = max(treadmilData.twoP);
    
    zeroAfterMax = find(treadmilData.twoP((maxTwoP_i + 1) : end) < maxTwoP_v) + maxTwoP_i;
    treadmilData(zeroAfterMax, :) = [];
    
    mkdir(fullfile(outputpath, 'Running'));
    mkdir(fullfile(outputpath, 'notRunning'));
    
    [filepath,~,~] = fileparts(videoPath);
    videoOutpath = fullfile(filepath, 'glm_vid');
    videoReader = VideoReader(videoPath);
    
%     splinesVal = {'0.05', '0.1', '0.25', '0.5', '0.75', '1', '2', '3', '6'};
    splinesVal = {'0.25', '0.5', '2'};
    
    for f_i = 1:length(splinesVal)
        splinesFile = sprintf('%s\\splines%s.csv', splinesLocation, splinesVal{f_i});
        
        splinesFuns{f_i}.func=xlsread(splinesFile);
        splinesFuns{f_i}.func = splinesFuns{f_i}.func(2:end, :);
        splinesFuns{f_i}.func = splinesFuns{f_i}.func(:, 2:end);
        splinesFuns{f_i}.delay = [0, 0.5, -0.5, 1, -1];
        splinesFuns{f_i}.name = replace(sprintf('s_%s', splinesVal{f_i}), '.', '');
    end
    
    behaveP = createPredictors(BehaveDataTreadmil, splinesFuns, ImageSamplingRate);
    
    indexLocationList = GLMRunning(doBehave, neuronNumberName, BehaveDataTreadmil, roiActivity, roiActivityNames, ImageSamplingRate, outputpath,videoOutpath, aftertonetime, treadmilData, videoSamplingRate, videoReader, behaveP); 
    GLMNotRunning(doBehave, neuronNumberName, BehaveDataTreadmil, roiActivity, roiActivityNames, ImageSamplingRate, outputpath,videoOutpath, aftertonetime, treadmilData, videoSamplingRate, videoReader, behaveP, indexLocationList);
end

function behaveP = createPredictors(BehaveDataTreadmil, splinesFuns, ImagingSamplingRate)
    speed_x2 = power(BehaveDataTreadmil.speed.*100, 2);
    speed_x3 = power(BehaveDataTreadmil.speed.*100, 3);
    accel_x2 = power(BehaveDataTreadmil.accel.*100, 2);
    accel_x3 = power(BehaveDataTreadmil.accel.*100, 3);
    binaryOnset = zeros(size(BehaveDataTreadmil.speed));
    binaryOffset = zeros(size(BehaveDataTreadmil.speed));   
    
    for i = 1:length(BehaveDataTreadmil.start)
        binaryOnset(max(BehaveDataTreadmil.start(i) - 5, 1):min(BehaveDataTreadmil.start(i) + 5, length(BehaveDataTreadmil.speed))) = 1;
        binaryOffset(max(BehaveDataTreadmil.stop(i) - 5, 1):min(BehaveDataTreadmil.stop(i) + 5, length(BehaveDataTreadmil.speed))) = 1;
    end
   
    behaveP.onset.reg.data = binaryOnset;
    behaveP.offset.reg.data = binaryOffset;
    behaveP.speed.reg.data = BehaveDataTreadmil.speed;
    behaveP.accel.reg.data = BehaveDataTreadmil.accel;
    behaveP.xLocation.reg.data = BehaveDataTreadmil.x_location;
    
    behaveP.speed_p2.reg.data = speed_x2;
    behaveP.speed_p3.reg.data = speed_x3;
    behaveP.accel_p2.reg.data = accel_x2;
    behaveP.accel_p3.reg.data = accel_x3;
    behaveP.posaccel.reg.data = BehaveDataTreadmil.posaccelBinary;
    behaveP.negaccel.reg.data = BehaveDataTreadmil.negaccelBinary;
    behaveP.walk.reg.data = BehaveDataTreadmil.walkBinary;
    behaveP.rest.reg.data = BehaveDataTreadmil.restBinary;    
    
    behaveP.speed.delay05.data = [zeros(ImagingSamplingRate*0.5, 1); BehaveDataTreadmil.speed(1:(end - ImagingSamplingRate*0.5))];
    behaveP.speed.early05.data = [BehaveDataTreadmil.speed((ImagingSamplingRate*0.5+1):end); zeros(ImagingSamplingRate*0.5, 1)];
    behaveP.accel.delay05.data = [zeros(ImagingSamplingRate*0.5, 1); BehaveDataTreadmil.accel(1:(end - ImagingSamplingRate*0.5))];
    behaveP.accel.early05.data = [BehaveDataTreadmil.accel((ImagingSamplingRate*0.5+1):end); zeros(ImagingSamplingRate*0.5, 1)];
    
    behaveP.speed_p2.delay05.data = [zeros( ImagingSamplingRate*0.5, 1); speed_x2(1:(end - ImagingSamplingRate*0.5))];
    behaveP.speed_p2.early05.data = [speed_x2((ImagingSamplingRate*0.5+1):end); zeros(ImagingSamplingRate*0.5, 1)];
    behaveP.accel_p2.delay05.data = [zeros( ImagingSamplingRate*0.5, 1); accel_x2(1:(end - ImagingSamplingRate*0.5))];
    behaveP.accel_p2.early05.data = [accel_x2((ImagingSamplingRate*0.5+1):end); zeros(ImagingSamplingRate*0.5, 1)];
    
    behaveP.speed_p3.delay05.data = [zeros( ImagingSamplingRate*0.5, 1); speed_x3(1:(end - ImagingSamplingRate*0.5))];
    behaveP.speed_p3.early05.data = [speed_x3((ImagingSamplingRate*0.5+1):end); zeros(ImagingSamplingRate*0.5, 1)];
    behaveP.accel_p3.delay05.data = [zeros( ImagingSamplingRate*0.5, 1); accel_x3(1:(end - ImagingSamplingRate*0.5))];
    behaveP.accel_p3.early05.data = [accel_x3((ImagingSamplingRate*0.5+1):end); zeros(ImagingSamplingRate*0.5, 1)];    
    
    behaveP.walk.delay05.data = [zeros( ImagingSamplingRate*0.5, 1); BehaveDataTreadmil.walkBinary(1:(end - ImagingSamplingRate*0.5))];
    behaveP.walk.early05.data = [BehaveDataTreadmil.walkBinary((ImagingSamplingRate*0.5+1):end); zeros( ImagingSamplingRate*0.5, 1)];
    behaveP.rest.delay05.data = [zeros( ImagingSamplingRate*0.5, 1); BehaveDataTreadmil.restBinary(1:(end - ImagingSamplingRate*0.5))];
    behaveP.rest.early05.data = [BehaveDataTreadmil.restBinary((ImagingSamplingRate*0.5+1):end); zeros( ImagingSamplingRate*0.5, 1)];
    
    
    fName = {'onset', 'offset', 'posaccel', 'negaccel'};
    
    for b_i = 1:length(fName)
        for sp_i = 1:length(splinesFuns)
            for k = 1:size(splinesFuns{sp_i}.func, 2)
                for j = 1:length(splinesFuns{sp_i}.delay)
                    delayF = ceil(abs(splinesFuns{sp_i}.delay(j)*ImagingSamplingRate));
                    if splinesFuns{sp_i}.delay(j) > 0
                        indicatorMatDelay = [zeros(delayF, 1); behaveP.(fName{b_i}).reg.data(1:(end-delayF))];
                    elseif splinesFuns{sp_i}.delay(j) < 0
                        indicatorMatDelay = [behaveP.(fName{b_i}).reg.data((delayF+1):end); zeros(delayF, 1)];
                    else
                        indicatorMatDelay = behaveP.(fName{b_i}).reg.data;
                    end

                    sp_del = splinesFuns{sp_i}.func(:, k);              
                    behaveP.(fName{b_i}).(sprintf('%s_%d_delay%d', splinesFuns{sp_i}.name, k, j)).data = filter(sp_del, 1, indicatorMatDelay);                   
                end
            end
        end
    end
end

function indexLocationList = GLMRunning(doBehave, neuronNumberName, BehaveDataTreadmil, roiActivity, roiActivityNames, ImageSamplingRate, outputpath,videoOutpath, aftertonetime, treadmilData, videoSamplingRate, videoReader, behaveP)   
    mkdir([videoOutpath, '\Running\']);
       
    TPAIndex = 1;
    indexLocationList = [];
    for onsetIndex = 1:length(BehaveDataTreadmil.start)
        indexLocationList  = [indexLocationList, (BehaveDataTreadmil.start(onsetIndex) - 2*ImageSamplingRate + 1):(BehaveDataTreadmil.stop(onsetIndex) + 2*ImageSamplingRate)];
    end
    
    indexLocationList = unique(indexLocationList);
    
    fNames = fieldnames(behaveP);
    
    for tr_i = 1:(aftertonetime+4)*ImageSamplingRate:length(indexLocationList)
        if (tr_i+(aftertonetime+4)*ImageSamplingRate-1) > length(indexLocationList)
            break;
        end
        
        currentLocation = indexLocationList(tr_i:(tr_i+(aftertonetime+4)*ImageSamplingRate-1));
        
        if any(currentLocation > size(roiActivity, 1)) || any(currentLocation > size(BehaveDataTreadmil.accel, 1))
            break;
        end
        
        createNewTPA(neuronNumberName,  roiActivity, roiActivityNames, ImageSamplingRate, aftertonetime, currentLocation, outputpath, 'Running', TPAIndex);            
        
        if doBehave
            createNewBDA(BehaveDataTreadmil, currentLocation, ImageSamplingRate, aftertonetime, aftertonetime+4, outputpath, 'Running', TPAIndex);

            for b_i = 1:length(fNames)
                peN = fieldnames(behaveP.(fNames{b_i}));
                for b2_i = 1:length(peN)
                    predictor.(fNames{b_i}).(peN{b2_i})(TPAIndex, :) = behaveP.(fNames{b_i}).(peN{b2_i}).data(currentLocation);
                end
            end

%             createVideo(treadmilData, currentLocation, videoSamplingRate, aftertonetime+4, videoReader, TPAIndex, 'Running', [videoOutpath, '\Running\']);
        end
        
        TPAIndex = TPAIndex + 1;
    end
            
    if doBehave
        save(fullfile(outputpath, 'Running', 'predictor_running'), 'predictor');
    end
end

function createNewBDA(BehaveDataTreadmil, indexLocationList, ImageSamplingRate, aftertonetime, trialSec, outputpath, Type, TPAIndex)
    strEvent = {};
    bdastrIndex = 1;
    staI = 1;
    stoI = 1;
    for eventI = 1:length(BehaveDataTreadmil.start)
        locI = find(indexLocationList == BehaveDataTreadmil.start(eventI));
        if ~isempty(locI)

            strEvent{bdastrIndex} = TPA_EventManager();
            strEvent{bdastrIndex}.Name = sprintf('%s:%02d', 'onset', staI);
            staI = staI + 1;
            strEvent{bdastrIndex}.tInd = [locI,locI +  5];
            strEvent{bdastrIndex}.Data = zeros(1, ImageSamplingRate*(aftertonetime+4));
            strEvent{bdastrIndex}.Data((locI):(locI +  5)) = 1;

            bdastrIndex = bdastrIndex + 1;
        end
    end
    
    for eventI = 1:length(BehaveDataTreadmil.stop)
        locI2 = find(indexLocationList == BehaveDataTreadmil.stop(eventI));
        if ~isempty(locI2)
            strEvent{bdastrIndex} = TPA_EventManager();
            strEvent{bdastrIndex}.Name = sprintf('%s:%02d', 'offset', stoI);
            stoI = stoI + 1;
            strEvent{bdastrIndex}.tInd = [locI2,locI2 +  5];
            strEvent{bdastrIndex}.Data = zeros(1, ImageSamplingRate*(aftertonetime+4));
            strEvent{bdastrIndex}.Data((locI2):(locI2 +  5)) = 1;

            bdastrIndex = bdastrIndex + 1;

        end
    end
  
    strEvent{bdastrIndex} = TPA_EventManager();
    strEvent{bdastrIndex}.Name = sprintf('posaccel');
    strEvent{bdastrIndex}.tInd = [1, ImageSamplingRate*(trialSec)];
    strEvent{bdastrIndex}.Data = BehaveDataTreadmil.posaccelBinary(indexLocationList);
    bdastrIndex = bdastrIndex + 1;

    strEvent{bdastrIndex} = TPA_EventManager();
    strEvent{bdastrIndex}.Name = sprintf('negaccel');
    strEvent{bdastrIndex}.tInd = [1, ImageSamplingRate*(trialSec)];
    strEvent{bdastrIndex}.Data = BehaveDataTreadmil.negaccelBinary(indexLocationList);
    bdastrIndex = bdastrIndex + 1;

    strEvent{bdastrIndex} = TPA_EventManager();
    strEvent{bdastrIndex}.Name = sprintf('speed');
    strEvent{bdastrIndex}.tInd = [1, ImageSamplingRate*(trialSec)];
    strEvent{bdastrIndex}.Data = BehaveDataTreadmil.speed(indexLocationList);
    bdastrIndex = bdastrIndex + 1;

    strEvent{bdastrIndex} = TPA_EventManager();
    strEvent{bdastrIndex}.Name = sprintf('accel');
    strEvent{bdastrIndex}.tInd = [1, ImageSamplingRate*(trialSec)];
    strEvent{bdastrIndex}.Data = BehaveDataTreadmil.accel(indexLocationList);      

    save(fullfile(outputpath, Type, sprintf('%s_%s_%03d', 'BDA', Type, TPAIndex)), 'strEvent');
end

function createNewTPA(neuronNumberName, roiActivity, roiActivityNames, ImageSamplingRate, aftertonetime, indexLocationList, outputpath, Type, TPAIndex)
    strROI = {};
    for roiIndex = 1:size(roiActivity, 2)
        strROI{roiIndex} = TPA_RoiManager();
        strROI{roiIndex}.Name = roiActivityNames{roiIndex};
        strROI{roiIndex}.Data = zeros(ImageSamplingRate*(aftertonetime+4), 2);
        strROI{roiIndex}.Data(:, 2) = roiActivity(indexLocationList, roiIndex);

        strROI{roiIndex}.xyInd = zeros(30, 2);
        strShift = zeros(ImageSamplingRate*(aftertonetime+4), 2);
    end

    save(fullfile(outputpath, Type, sprintf('%s_%s_%s_%03d', 'TPA',neuronNumberName, Type, TPAIndex)), 'strROI', 'strShift');
end

function GLMNotRunning(doBehave, neuronNumberName, BehaveDataTreadmil, roiActivity, roiActivityNames, ImageSamplingRate, outputpath,videoOutpath, aftertonetime, treadmilData, videoSamplingRate, videoReader, behaveP, indexLocationListRunning)   
    mkdir([videoOutpath, '\notRunning\']);
    TPAIndex = 1;
    
    indexLocationList = 1:length(BehaveDataTreadmil.speed);
    
    indexLocationListRunning(indexLocationListRunning > length(BehaveDataTreadmil.speed)) = [];
    indexLocationList(indexLocationListRunning) = [];
    
    fNames = fieldnames(behaveP);
    
    for tr_i = 1:(aftertonetime+4)*ImageSamplingRate:length(indexLocationList)
        if (tr_i+(aftertonetime+4)*ImageSamplingRate-1) > length(indexLocationList)
            break;
        end
        
        currentLocation = indexLocationList(tr_i:(tr_i + (aftertonetime+4)*ImageSamplingRate - 1));
        
        if any(currentLocation > size(roiActivity, 1))
            break;
        end
        
        
        createNewTPA(neuronNumberName, roiActivity, roiActivityNames, ImageSamplingRate, aftertonetime, currentLocation, outputpath, 'notRunning', TPAIndex);            
        
        if doBehave
            createNewBDA(BehaveDataTreadmil, currentLocation, ImageSamplingRate, aftertonetime, aftertonetime+4, outputpath, 'notRunning', TPAIndex);           

            for b_i = 1:length(fNames)
                peN = fieldnames(behaveP.(fNames{b_i}));
                for b2_i = 1:length(peN)
                    predictor.(fNames{b_i}).(peN{b2_i})(TPAIndex, :) = behaveP.(fNames{b_i}).(peN{b2_i}).data(currentLocation);
                end
            end

%             createVideo(treadmilData, currentLocation, videoSamplingRate, aftertonetime+4, videoReader, TPAIndex, 'notRunning',  [videoOutpath, '\notRunning\']);
        end
        TPAIndex = TPAIndex + 1;
    end
    
    if doBehave
        save(fullfile(outputpath, 'notRunning', 'predictor_notRunning'), 'predictor');
    end
end

function createVideo(treadmilData, indexLocationList, videoSamplingRate, trialInsec, videoReader, TPA_index, TPA_type, outputpath)
    firstFrameOffset = treadmilData.twoP(1);
    treadmilData.twoP = treadmilData.twoP - firstFrameOffset;
    
    firstF = 0;
    for i = 1:length(indexLocationList)
        firstFrameFromVideo = find(treadmilData.twoP == indexLocationList(i)-1, 1, 'first');
        lastFrameFromVideo = find(treadmilData.twoP == indexLocationList(i)-1, 1, 'last');
   
        if isempty(firstFrameFromVideo)
            continue;
        end
        
        if ~firstF
            frame = read(videoReader,[firstFrameFromVideo,lastFrameFromVideo]);
            firstF = 1;
        else
            frame(:,:,:, end+1:(end+lastFrameFromVideo-firstFrameFromVideo+1)) = read(videoReader,[firstFrameFromVideo,lastFrameFromVideo]);     
        end
    end
     
    v = VideoWriter(fullfile(outputpath, sprintf('Movie_%03d_%s', TPA_index, TPA_type)));
    
    totalFrames = trialInsec*videoSamplingRate;
    for i = size(frame, 4)+1:totalFrames
        frame(:, :, :, i) = frame(:, :, :, end);
    end
  
    
    open(v);
    writeVideo(v,frame)
    close(v)
end