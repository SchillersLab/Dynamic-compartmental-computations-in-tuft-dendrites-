function [SpikeTrainStart, SpikeTrainEnd, SpikeTrainPks,SpikeTrainH, eventDetector_activity, comb_activity] = calcRoiEventDetectorByMLSpike_V3_Sim(dataCurROI, ImageSampleTime, frameNum, aV, outputpath, activityIndex, clusterCount, roiName, sigmaChangeValue, roiCount, runMLS, thresholdForGn)       
    traceSig            = dataCurROI;
        
    sampFreq        = 1/ImageSampleTime;
    
    SpikeTrainStart  = [];
    SpikeTrainEnd  = [];
    SpikeTrainPks  = [];
    SpikeTrainH = [];
    max_index = [];
    
    [pksV,locpks, wpks] = findpeaks(traceSig,'MinPeakHeight', thresholdForGn);
    
    slopTh = diff(traceSig);
    slopTh = [slopTh(1); slopTh];
    for i = 1:length(pksV)
        startLocation = find(slopTh(1:(locpks(i)-1)) <= 0, 1, 'last');
        endLocation = find(slopTh((locpks(i)+1):end) >= 0, 1, 'first');
        
        if isempty(endLocation)
            SpikeTrainEnd(i) = length(traceSig);
        else
            SpikeTrainEnd(i) = endLocation + locpks(i);
        end
        if isempty(startLocation)
            startLocation = 1;
        else
            SpikeTrainStart(i) =  startLocation;
        end
        SpikeTrainPks(i) = locpks(i);
        SpikeTrainH(i) = pksV(i);
    end
    
                 
    for i = 1:length(SpikeTrainStart)
        tr_indexStart = floor(SpikeTrainStart(i) ./ frameNum) + 1;
        tr_indexPKS = floor(SpikeTrainPks(i) ./ frameNum) + 1;
       
        nextLoc = (tr_indexPKS - 1) * frameNum + 1;
        if tr_indexStart ~= tr_indexPKS & SpikeTrainPks(i) >= nextLoc
             SpikeTrainStart(i) = nextLoc;
        end
    end
    
    
    f = figure;
    hold on;
    
    sb1 = subplot(8, 1, 1:6);
    hold on;
    title(roiName);
    
    plot(traceSig)
    
    plot(SpikeTrainPks, traceSig(SpikeTrainPks), '*r');
    plot(SpikeTrainStart, traceSig(SpikeTrainStart), '*b');
    
    legend('Activity', 'Fit', 'Peaks', 'StartEvent', 'EndEvents');
    
    plot(SpikeTrainPks, SpikeTrainH, 'o');
    
    xlim([1, size(traceSig, 1)]);
%     ylim([-1, 5]);
    
    sb2 = subplot(8, 1, 8:8);
    imagesc(traceSig');
    colormap(jet);
    caxis(sb1.YLim);
       
    linkaxes([sb1, sb2], 'x');
    
    mysave(f, [outputpath, '\activity_averagePksHistByMLSpike_', num2str(activityIndex)]);
    
    close(f);
    
    comb_activity = traceSig;
    eventDetector_activity = zeros(length(traceSig), 1);
    for i = 1:length(SpikeTrainPks)
        
        eventDetector_activity((SpikeTrainStart(i)):(SpikeTrainPks(i))) = ones(length((SpikeTrainStart(i)):SpikeTrainPks(i)),1);
    end
end