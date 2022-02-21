function saveDataFotGLMTreadMil_Type2(BehaveDataTreadmil, roiActivityLength, ImageSamplingRate, outputpath, trialTime,treadMilTxt, videoPath)
    treadmilData = readtable(treadMilTxt);   
    
    mkdir(fullfile(outputpath, 'BDA_new'));
    
    mkdir(fullfile(fileparts(videoPath), 'NewVideos'));
    
    videoReader = VideoReader(videoPath);
    
    TPAIndex = 1;
    startTwoP = treadmilData.twoP(1);
    
    for imagingIndex = 1:trialTime*ImageSamplingRate:roiActivityLength
        startIndexImaging(TPAIndex) = imagingIndex;
        endIndexImaging(TPAIndex) = imagingIndex + trialTime*ImageSamplingRate - 1;
        
        startIndexBehave(TPAIndex) = find(treadmilData.twoP == (startIndexImaging(TPAIndex) - 1 + startTwoP), 1, 'first');
        endIndexBehave(TPAIndex) = find(treadmilData.twoP == (endIndexImaging(TPAIndex) - 1 + startTwoP), 1, 'last');
%         
        TPAIndex = TPAIndex + 1;
    end
    
    trialSec = trialTime;
    
    for tr_i = 1:length(startIndexImaging)
        indexLocationList = startIndexImaging(tr_i) : endIndexImaging(tr_i);
     
        strEvent = {};
        bdastrIndex = 1;
        staI = 1;
        for eventI = 1:length(BehaveDataTreadmil.start)
            locI = find(indexLocationList == BehaveDataTreadmil.start(eventI), 1);
            if ~isempty(locI)
                strEvent{bdastrIndex} = TPA_EventManager();
                strEvent{bdastrIndex}.Name = sprintf('%s:%02d', 'onset', staI);
                staI = staI + 1;
                strEvent{bdastrIndex}.tInd = [max(locI - 5, 1),min(locI +  5, ImageSamplingRate*(trialSec))];
                strEvent{bdastrIndex}.Data = zeros(1, ImageSamplingRate*(trialSec));
                strEvent{bdastrIndex}.Data(max(locI - 5, 1):min(locI +  5, ImageSamplingRate*(trialSec))) = 1;
               
                bdastrIndex = bdastrIndex + 1;
            end
        end
        
        stoI = 1;
        for eventI = 1:length(BehaveDataTreadmil.stop)
            locI = find(indexLocationList == BehaveDataTreadmil.stop(eventI), 1);
            if ~isempty(locI)
                strEvent{bdastrIndex} = TPA_EventManager();
                strEvent{bdastrIndex}.Name = sprintf('%s:%02d', 'offset', stoI);
                stoI = stoI + 1;
                strEvent{bdastrIndex}.tInd = [max(locI - 5, 1),min(locI +  5, ImageSamplingRate*(trialSec))];
                strEvent{bdastrIndex}.Data = zeros(1, ImageSamplingRate*(trialSec));
                strEvent{bdastrIndex}.Data(max(locI - 5, 1):min(locI +  5, ImageSamplingRate*(trialSec))) = 1;
                
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
        bdastrIndex = bdastrIndex + 1;       
        
        save(fullfile(outputpath, 'BDA_new', sprintf('%s_%03d', 'BDA_Reg', tr_i)), 'strEvent');
        
        predictor.speed(tr_i, :) = BehaveDataTreadmil.speed(indexLocationList);
        predictor.accel(tr_i, :) = BehaveDataTreadmil.accel(indexLocationList);
        predictor.walk(tr_i, :) = BehaveDataTreadmil.walkBinary(indexLocationList);
        predictor.rest(tr_i, :) = BehaveDataTreadmil.restBinary(indexLocationList);
        predictor.posaccel(tr_i, :) = BehaveDataTreadmil.posaccelBinary(indexLocationList);
        predictor.negaccel(tr_i, :) = BehaveDataTreadmil.negaccelBinary(indexLocationList);
        predictor.x_location(tr_i, :) = BehaveDataTreadmil.x_location(indexLocationList);
        createVideo(startIndexBehave(tr_i), endIndexBehave(tr_i), videoReader, tr_i, [fileparts(videoPath), '\NewVideos\'], BehaveDataTreadmil.BehaveSamplingRate*trialTime);
    end
    
    save(fullfile(outputpath, 'BDA_new', 'predictor_reg'), 'predictor'); 
    
end

function createVideo(startIndexBehave, endIndexBehave, videoReader, tr_i, outputpath, totalFrames)
    firstFrameFromVideo = startIndexBehave;
    lastFrameFromVideo = endIndexBehave;
    
    v = VideoWriter(fullfile(outputpath, sprintf('Movie_%03d_%s', tr_i, 'Reg')));
    frame = read(videoReader,[firstFrameFromVideo,lastFrameFromVideo]); 
    
    for i = size(frame, 4)+1:totalFrames
        frame(:, :, :, i) = frame(:, :, :, end);
    end
    
    open(v);
    writeVideo(v,frame)
    close(v);
end