function [SpikeTrainStart, SpikeTrainEnd, SpikeTrainPks,SpikeTrainH, eventDetector_activity, comb_activity] = calcRoiEventDetectorByMLSpike_V3(dataCurROI, ImageSampleTime, frameNum, aV, outputpath, activityIndex, clusterCount, roiName, sigmaChangeValue, roiCount, runMLS, thresholdForGn)       
    traceSig            = dataCurROI;
        
    sampFreq        = 1/ImageSampleTime;
    
    SpikeTrainStart  = [];
    SpikeTrainEnd  = [];
    SpikeTrainPks  = [];
    SpikeTrainH = [];
    max_index = [];
    
    
    if runMLS
        par = tps_mlspikes('par');
        par.dt = ImageSampleTime;
        par.a = aV;
        if sigmaChangeValue ~= 0
            par.finetune.sigma = sigmaChangeValue;
        end
        par.drift.parameter = .015;
        par.dographsummary = false;

        [spk, fit, drift, parest] = spk_est(traceSig,par);

    % -------------------------------------------------------------------------------    
        countVec = fn_timevector(spk, ImageSampleTime, 'count');
    else
        spkmin = thresholdForGn*GetSn(traceSig);
        model_ar = 'ar2';
        fr = ImageSampleTime;
        decay_time = 0.5;
        lam = choose_lambda(exp(-1/(fr*decay_time)),GetSn(traceSig),0.99);
        [cc,spk,opts_oasis] = deconvolveCa(traceSig,model_ar,'method','thresholded','optimize_pars',true,'maxIter',20,...
                                    'window',150,'lambda',lam,'smin',spkmin);
        fit = cc;                        
        countVec = spk ~= 0;
    end
    
    countVec((end + 1) : length(traceSig)) = 0;
    
    isRunNeeded = 1;
    
    while isRunNeeded
        textSearch = num2str(countVec');
        textSearch = replace(textSearch, ' ', '');
        textSR_1 = regexp(num2str(textSearch),  '[1-9][0][1-9]');
        textSR_2 = regexp(num2str(textSearch),  '[1-9][0][0][1-9]');
        textSR_3 = regexp(num2str(textSearch),  '[1-9][0][0][0][1-9]');
        
        if isempty(textSR_1) && isempty(textSR_2) && isempty(textSR_3)
            isRunNeeded = 0;
        else 
            if ~isempty(textSR_1)
                countVec(textSR_1 + 1) = 1;
            end
            
            if ~isempty(textSR_2)
                countVec(textSR_2 + 1) = 1;
                countVec(textSR_2 + 2) = 1;
            end
            if ~isempty(textSR_3)
                countVec(textSR_3 + 1) = 1;
                countVec(textSR_3 + 2) = 1;
                countVec(textSR_3 + 3) = 1;
            end
            
        end
    end
    
    
    countVecFix = find(countVec(2:end) ~= 0 & countVec(1:end-1) == 0);
    
    v_calc = (fit(2:(end)) - fit(1:(end-1)));
    v_calc = repelem(v_calc,1,1);
    v_calc = [v_calc(1) ;v_calc];
   
    filtDur         = 0.6 * sampFreq;      % filter duration in sec
    filtLenH        = ceil(filtDur/2);
    filtLen         = filtLenH*2;

    filtSmooth          = hamming(filtLen);
    filtSmooth          = filtSmooth./sum(filtSmooth);

    filterRes = filtfilt(filtSmooth,1,traceSig);   
        
    
    v_calc_org = (filterRes(2:(end)) - filterRes(1:(end-1)));
    v_calc_org = repelem(v_calc_org,1,1);
    v_calc_org = [v_calc_org(1) ;v_calc_org];
   
    v_slop = v_calc_org > 0;
    res = [0; find(v_slop(1:end-1) ~= v_slop(2:end)); length(v_calc)];
    
    for k = 1:length(res)-1
        if v_slop(res(k) + 1) == 0 & (res(k + 1) - res(k) <= 1)
            slop_value = 1;
            if (k-1>0)
                slop_value = v_slop(res(k-1) + 1) ;    
            end

            v_slop((res(k) + 1):res(k+1)) = slop_value;
        end
    end
    
    for k = 1:length(res)-1
        if v_slop(res(k) + 1) == 1 & (res(k + 1) - res(k) <= 1)
            slop_value = 0;
            if (k-1>0)
                slop_value = v_slop(res(k-1) + 1) ;    
            end

            v_slop((res(k) + 1):res(k+1)) = slop_value;
        end
    end


    
    spk_index = 1;
    
    
    for i = 1:length(countVecFix)
        start_pos = countVecFix(i);
        
        nextPos = length(countVec);
        if i < length(countVecFix)
            nextPos = countVecFix(i + 1);
        end
        
        if (nextPos - start_pos <= 3)
            continue;
        end
        
        
        endPos = find(countVec((start_pos+1):nextPos) == 0, 1);
        endPos = endPos + start_pos;
        
        
        if isempty(endPos) || (endPos - start_pos == 3 & all(fit((start_pos + 1) : (endPos - 1)) == fit(start_pos + 1)))
            continue
        end
        
        peaksSub = find(v_slop(start_pos:(endPos-1)) == 1 & v_slop((start_pos+1):endPos) == 0);
        peaksSub = peaksSub + start_pos - 1;
        
        if length(peaksSub) >= 1
            
            if all(filterRes(peaksSub(1):peaksSub(end)) > min(fit(peaksSub(1):peaksSub(end))))
                start_pos_sub = start_pos;
                for j = 1:(length(peaksSub)-1)
                    maxPos = peaksSub(j);
                    naxPosSub = peaksSub(j + 1);
                    endPos2 = find(v_slop(maxPos:naxPosSub) == 0, 1, 'last');
                    endPos2 = endPos2 + maxPos - 1;

                    SpikeTrainStart(spk_index) = start_pos_sub;
                    SpikeTrainEnd(spk_index) = endPos2;
                    SpikeTrainPks(spk_index) = maxPos;
                    
                    [SpikeTrainH(spk_index), max_index(spk_index)] = max(traceSig(SpikeTrainStart(spk_index):SpikeTrainEnd(spk_index))); 
                    max_index(spk_index) = max_index(spk_index) + SpikeTrainStart(spk_index) - 1;       
                    spk_index = spk_index + 1; 
                    start_pos_sub = endPos2;               
                end

                maxPos = peaksSub(end);        
                start_pos = start_pos_sub;
            else
                [~, maxPos] = max(fit(start_pos:endPos));
                maxPos = maxPos + start_pos - 1;
            end
                        
        else
            [~, maxPos] = max(fit(start_pos:endPos));
            maxPos = maxPos + start_pos - 1;
        end
        
        endPos2 = find(countVec(start_pos:nextPos) == 0, 1, 'last');
        endPos2 = endPos2 + start_pos - 1;
        
        afterS = find(fit(endPos:endPos2) <= 0.9 * fit(maxPos(1)), 1);
        afterS = afterS + endPos - 1;
        end_find_res = find(abs(v_calc(afterS:nextPos)) < 0.0001, 1);
        
        if isempty(end_find_res)
            end_find_res = nextPos; 
        else
            end_find_res = end_find_res + afterS - 1;
        end
%         
        
        SpikeTrainStart(spk_index) = start_pos;
        SpikeTrainEnd(spk_index) = end_find_res;
        SpikeTrainPks(spk_index) = maxPos(1);
        
        [SpikeTrainH(spk_index), max_index(spk_index)] = max(traceSig(SpikeTrainStart(spk_index):SpikeTrainEnd(spk_index))); 
        max_index(spk_index) = max_index(spk_index) + SpikeTrainStart(spk_index) - 1;       
        spk_index = spk_index + 1;     
    end

    for i = 1:length(SpikeTrainStart)
        tr_indexStart = floor(SpikeTrainStart(i) ./ frameNum) + 1;
        tr_indexPKS = floor(SpikeTrainPks(i) ./ frameNum) + 1;
       
        if tr_indexStart ~= tr_indexPKS
             SpikeTrainStart(i) = (tr_indexPKS - 1) * frameNum + 1;
        end
    end
    
    
    f = figure;
    hold on;
    
    sb1 = subplot(8, 1, 1:6);
    hold on;
    title(roiName);
    
    plot(traceSig)
    plot(fit)

    plot(SpikeTrainPks, fit(SpikeTrainPks), '*r');
    plot(SpikeTrainStart, fit(SpikeTrainStart), '*b');
    plot(SpikeTrainEnd, fit(SpikeTrainEnd), '*g');
 
    legend('Activity', 'Fit', 'Peaks', 'StartEvent', 'EndEvents');
    
    plot(max_index, SpikeTrainH, 'o');
    
    xlim([1, size(traceSig, 1)]);
%     ylim([-1, 5]);
    
    sb2 = subplot(8, 1, 8:8);
    imagesc(traceSig');
    colormap(jet);
    caxis(sb1.YLim);
       
    linkaxes([sb1, sb2], 'x');
    
    mysave(f, [outputpath, '\activity_averagePksHistByMLSpike_', num2str(activityIndex)]);
    
    close(f);
    
    comb_activity = fit;
    eventDetector_activity = zeros(length(fit), 1);
    for i = 1:length(SpikeTrainPks)
        comb_activity(SpikeTrainStart(i):SpikeTrainEnd(i)) = traceSig(SpikeTrainStart(i):SpikeTrainEnd(i));
        
%         eventDetector_activity((SpikeTrainStart(i)):(SpikeTrainPks(i))) = ones(length((SpikeTrainStart(i)):SpikeTrainPks(i)),1) .* roiCount;
        
        eventDetector_activity((SpikeTrainStart(i)):(SpikeTrainPks(i))) = ones(length((SpikeTrainStart(i)):SpikeTrainPks(i)),1);
    end
end