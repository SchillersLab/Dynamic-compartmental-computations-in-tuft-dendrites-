function saveDataFotGLMTreadMil(BehaveDataTreadmil, roiActivity, roiActivityNames, ImageSamplingRate, tpaFolder, aftertonetime, treadMilTxt, videoPath)
    treadmilData = readtable(treadMilTxt);   
    
    [filepath,name,~] = fileparts(tpaFolder);
    
    outputpath = fullfile(filepath, [name, '_glm']);
    
    videoTimeInsec = ((treadmilData.time(end) - treadmilData.time(1)) ./  1000);
    videoSamplingRate = ceil(size(treadmilData, 1) ./ videoTimeInsec);
    
    mkdir(fullfile(outputpath, 'Onset'));
    mkdir(fullfile(outputpath, 'Offset'));
    
    [filepath,~,~] = fileparts(videoPath);
    videoOutpath = fullfile(filepath, 'glm_vid');
    videoReader = VideoReader(videoPath);
    
    GLMOnset(BehaveDataTreadmil, roiActivity, roiActivityNames, ImageSamplingRate, outputpath,videoOutpath, aftertonetime, treadmilData, videoSamplingRate, videoReader); 
    GLMOffset(BehaveDataTreadmil, roiActivity, roiActivityNames, ImageSamplingRate, outputpath,videoOutpath, aftertonetime, treadmilData, videoSamplingRate, videoReader);
end

function GLMOnset(BehaveDataTreadmil, roiActivity, roiActivityNames, ImageSamplingRate, outputpath,videoOutpath, aftertonetime, treadmilData, videoSamplingRate, videoReader)   
    mkdir([videoOutpath, '\Onset\']);
       
    TPAIndex = 1;
    for onsetIndex = 1:length(BehaveDataTreadmil.start)
        if (onsetIndex == 1 && BehaveDataTreadmil.start(onsetIndex) >= ImageSamplingRate*4) || ...
            (onsetIndex ~= 1 && ...
                (BehaveDataTreadmil.start(onsetIndex) - BehaveDataTreadmil.start(onsetIndex - 1) >=  ImageSamplingRate*4))
            
            indexLocationList  = (BehaveDataTreadmil.start(onsetIndex) - 4*ImageSamplingRate + 1):(BehaveDataTreadmil.start(onsetIndex) + aftertonetime*ImageSamplingRate);
              
            if indexLocationList(end) > size(roiActivity, 1)
                continue;
            end
            
            createNewTPA(roiActivity, roiActivityNames, ImageSamplingRate, aftertonetime, indexLocationList, outputpath, 'Onset', TPAIndex);            
            createNewBDA(BehaveDataTreadmil, indexLocationList, ImageSamplingRate, aftertonetime, aftertonetime+4, outputpath, 'Onset', TPAIndex);
            
            predictor.speed(TPAIndex, :) = BehaveDataTreadmil.speed(indexLocationList);
            predictor.accel(TPAIndex, :) = BehaveDataTreadmil.accel(indexLocationList);
            predictor.walk(TPAIndex, :) = BehaveDataTreadmil.walkBinary(indexLocationList);
            predictor.rest(TPAIndex, :) = BehaveDataTreadmil.restBinary(indexLocationList);
            predictor.posaccel(TPAIndex, :) = BehaveDataTreadmil.posaccelBinary(indexLocationList);
            predictor.negaccel(TPAIndex, :) = BehaveDataTreadmil.negaccelBinary(indexLocationList);
            predictor.x_location(TPAIndex, :) = BehaveDataTreadmil.x_location(indexLocationList);
            
            createVideo(treadmilData, indexLocationList, videoSamplingRate, aftertonetime+4, videoReader, TPAIndex, 'onset', [videoOutpath, '\Onset\']);
            
            TPAIndex = TPAIndex + 1;
        end
    end
    
    save(fullfile(outputpath, 'Onset', 'predictor_onset'), 'predictor');            
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

function createNewTPA(roiActivity, roiActivityNames, ImageSamplingRate, aftertonetime, indexLocationList, outputpath, Type, TPAIndex)
    strROI = {};
    for roiIndex = 1:size(roiActivity, 2)
        strROI{roiIndex} = TPA_RoiManager();
        strROI{roiIndex}.Name = roiActivityNames{roiIndex};
        strROI{roiIndex}.Data = zeros(ImageSamplingRate*(aftertonetime+4), 2);
        strROI{roiIndex}.Data(:, 2) = roiActivity(indexLocationList, roiIndex);

        strROI{roiIndex}.xyInd = zeros(30, 2);
        strShift = zeros(ImageSamplingRate*(aftertonetime+4), 2);
    end

    save(fullfile(outputpath, Type, sprintf('%s_%s_%03d', 'TPA', Type, TPAIndex)), 'strROI', 'strShift');
end

function GLMOffset(BehaveDataTreadmil, roiActivity, roiActivityNames, ImageSamplingRate, outputpath,videoOutpath, aftertonetime, treadmilData, videoSamplingRate, videoReader)   
    mkdir([videoOutpath, '\Offset\']);
    TPAIndex = 1;
    for offsetIndex = 1:length(BehaveDataTreadmil.stop)
        if (offsetIndex == 1 && BehaveDataTreadmil.stop(offsetIndex) >= ImageSamplingRate*4) || ...
            (offsetIndex ~= 1 && ...
                (BehaveDataTreadmil.stop(offsetIndex) - BehaveDataTreadmil.stop(offsetIndex - 1) >=  ImageSamplingRate*4))
            
            indexLocationList  = (BehaveDataTreadmil.stop(offsetIndex) - 4*ImageSamplingRate + 1):(BehaveDataTreadmil.stop(offsetIndex) + aftertonetime*ImageSamplingRate);
               
            if indexLocationList(end) > size(roiActivity, 1)
                continue;
            end
            
            createNewTPA(roiActivity, roiActivityNames, ImageSamplingRate, aftertonetime, indexLocationList, outputpath, 'Offset', TPAIndex);            
            createNewBDA(BehaveDataTreadmil, indexLocationList, ImageSamplingRate, aftertonetime, aftertonetime+4, outputpath, 'Offset', TPAIndex);           
            
            predictor.speed(TPAIndex, :) = BehaveDataTreadmil.speed(indexLocationList);
            predictor.accel(TPAIndex, :) = BehaveDataTreadmil.accel(indexLocationList);
            predictor.walk(TPAIndex, :) = BehaveDataTreadmil.walkBinary(indexLocationList);
            predictor.rest(TPAIndex, :) = BehaveDataTreadmil.restBinary(indexLocationList);
            predictor.x_location(TPAIndex, :) = BehaveDataTreadmil.x_location(indexLocationList);
            predictor.posaccel(TPAIndex, :) = BehaveDataTreadmil.posaccelBinary(indexLocationList);
            predictor.negaccel(TPAIndex, :) = BehaveDataTreadmil.negaccelBinary(indexLocationList);
            
            createVideo(treadmilData, indexLocationList, videoSamplingRate, aftertonetime+4, videoReader, TPAIndex, 'offset',  [videoOutpath, '\Offset\']);
            
            TPAIndex = TPAIndex + 1;
        end
    end
    
    save(fullfile(outputpath, 'Offset', 'predictor_offset'), 'predictor');            
end

function createVideo(treadmilData, indexLocationList, videoSamplingRate, trialInsec, videoReader, TPA_index, TPA_type, outputpath)
    firstFrameOffset = treadmilData.twoP(1);
    treadmilData.twoP = treadmilData.twoP - firstFrameOffset;
    firstFrameFromVideo = find(treadmilData.twoP == indexLocationList(1)-1, 1, 'first');
    lastFrameFromVideo = find(treadmilData.twoP == indexLocationList(end)-1, 1, 'last');
    
    v = VideoWriter(fullfile(outputpath, sprintf('Movie_%03d_%s', TPA_index, TPA_type)));
    frame = read(videoReader,[firstFrameFromVideo,lastFrameFromVideo]); 
    
    totalFrames = trialInsec*videoSamplingRate;
    for i = size(frame, 4)+1:totalFrames
        frame(:, :, :, i) = frame(:, :, :, end);
    end
  
    
    open(v);
    writeVideo(v,frame)
    close(v)
end