function [SpikeTrainStart, SpikeTrainEnd, SpikeTrainPks, SpikeTrainH, SpikeTrainClusterSec, comb_activity] = calcRoiEventDetectorByMLSpike(dataCurROI, ImageSampleTime, frameNum, aV, outputpath, activityIndex, clusterCount, roiName, sigmaChangeValue)       
    traceSig            = dataCurROI;
        
    sampFreq        = 1/ImageSampleTime;
    
    SpikeTrainStart  = [];
    SpikeTrainEnd  = [];
    SpikeTrainPks  = [];
    SpikeTrainH = [];
    max_index = [];
    SpikeTrainCluster = [];
%     
    % parameters
    par = tps_mlspikes('par');
    par.dt = ImageSampleTime;
%     par.algo.estimate = 'proba';
    % (use autocalibrated parameters)
    par.a = aV;
%     
%     par.tau = 0.6;

    if sigmaChangeValue ~= 0
        par.finetune.sigma = sigmaChangeValue;
    end
%     

% (the OGB saturation and drift parameters are fixed)
%     par.saturation = 0.1;
    par.drift.parameter = .015;
    % (do not display graph summary)
%     par.dographsummary = false;

    [spk, fit, drift, parest] = spk_est(traceSig,par);
    
% -------------------------------------------------------------------------------    
    countVec = fn_timevector(spk, ImageSampleTime, 'count');
    
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
        
        
        if endPos - start_pos == 3 & all(fit((start_pos + 1) : (endPos - 1)) == fit(start_pos + 1))
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
        
        afterS = find(fit(endPos:endPos2) <= 0.9 * fit(maxPos), 1);
        afterS = afterS + endPos - 1;
        end_find_res = find(abs(v_calc(afterS:nextPos)) < 0.0001, 1);
        
        if isempty(end_find_res)
            end_find_res = nextPos; 
        else
            end_find_res = end_find_res + afterS - 1;
        end
%         
%         if end_find_res - start_pos <= 4
%             continue;
%         end
        
        SpikeTrainStart(spk_index) = start_pos;
        SpikeTrainEnd(spk_index) = end_find_res;
        SpikeTrainPks(spk_index) = maxPos;
        
        [SpikeTrainH(spk_index), max_index(spk_index)] = max(traceSig(SpikeTrainStart(spk_index):SpikeTrainEnd(spk_index))); 
        max_index(spk_index) = max_index(spk_index) + SpikeTrainStart(spk_index) - 1;       
        spk_index = spk_index + 1;     
    end
    
% ---------------------------------------------------------------------------------    
% 
%     
%     spk = round(spk * sampFreq);  
%     
%     if isempty(spk)
%         return;
%     end
%     
%     spk = unique(spk);
%     spk = [-1, spk];  
%     spk_fin = spk(find(spk(2:end) - spk(1:(end - 1)) > 4) + 1);
%     
%     v_calc = (fit(2:(end)) - fit(1:(end-1)));
%     v_calc = repelem(v_calc,1,1);
%     v_calc = [v_calc(1) ;v_calc];
%    
%     v_slop = v_calc > 0;
% 
%     res = [0; find(v_slop(1:end-1) ~= v_slop(2:end)); length(v_calc)];
% 
%     for k = 1:length(res)-1
%         if v_slop(res(k) + 1) == 0 & (res(k + 1) - res(k) <= 1)
%             slop_value = 1;
%             if (k-1>0)
%                 slop_value = v_slop(res(k-1) + 1) ;    
%             end
% 
%             v_slop((res(k) + 1):res(k+1)) = slop_value;
%             
%             if slop_value == 1
%                 v_calc((res(k) + 1):res(k+1)) = 0.1;
%             else
%                 v_calc((res(k) + 1):res(k+1)) = -0.1;
%             end
%         end
%     end
% 
%     
%     spk_index = 1;
%     
%     for i = 1:length(spk_fin)
%         start_pos = spk_fin(i);
%         
%         nextPos = length(fit);
%         if i < length(spk_fin)
%             nextPos = spk_fin(i + 1);
%         end
%         
%         maxRes = find(fit((start_pos+1):nextPos) <= fit(start_pos), 1);
%         if isempty(maxRes)
%             maxRes = nextPos; 
%         else
%             maxRes = maxRes(1) + start_pos;
%         end
%         
%         maxRes2 = maxRes;
%         maxRes3 = find(fit(start_pos:maxRes2) == max(fit(start_pos:maxRes2)), 1);
% 
%         pks_Pos_n = maxRes3(1) + start_pos - 1;
%                                
%         afterS = find(fit(pks_Pos_n:nextPos) <= 0.9 * fit(pks_Pos_n), 1);
%         afterS = afterS + pks_Pos_n - 1;
%         end_find_res = find(abs(v_calc(afterS:nextPos)) < 0.0001, 1);
%         
%         if isempty(end_find_res)
%             end_find_res = nextPos; 
%         else
%             end_find_res = end_find_res + afterS - 1;
%         end
%         
% %          peaksSub = find(v_slop(start_pos:(end_find_res-1)) == 1 & v_slop((start_pos+1):end_find_res) == 0);
% %         
% %         if length(peaksSub) >= 1
% %             
% %             start_pos_sub = start_pos;
% %             for j = 1:(length(peaksSub)-1)
% %                 maxPos = peaksSub(j) + start_pos - 1;
% %                 naxPosSub = peaksSub(j + 1) + start_pos - 1;
% %                 endPos2 = find(v_slop(maxPos:naxPosSub) == 0, 1, 'last');
% %                 endPos2 = endPos2 + maxPos - 1;
% %                 
% %                 SpikeTrainStart(spk_index) = start_pos_sub;
% %                 SpikeTrainEnd(spk_index) = endPos2;
% %                 SpikeTrainPks(spk_index) = maxPos;
% % 
% %                 [SpikeTrainH(spk_index), max_index(spk_index)] = max(traceSig(SpikeTrainStart(spk_index):SpikeTrainEnd(spk_index))); 
% %                 max_index(spk_index) = max_index(spk_index) + SpikeTrainStart(spk_index) - 1;       
% %                 spk_index = spk_index + 1; 
% %                 start_pos_sub = endPos2;               
% %             end
% %             
% %             maxPos = peaksSub(end) + start_pos - 1;        
% %             start_pos = start_pos_sub; 
% %             pks_Pos_n = maxPos;
% %         end
%        
%         
%         SpikeTrainStart(spk_index) = start_pos;
%         SpikeTrainPks(spk_index) = pks_Pos_n;
%         SpikeTrainEnd(spk_index) = end_find_res;
%         
%         [SpikeTrainH(spk_index), max_index(spk_index)] = max(traceSig(SpikeTrainStart(spk_index):SpikeTrainEnd(spk_index))); 
%          max_index(spk_index) = max_index(spk_index) + SpikeTrainStart(spk_index) - 1;       
%         spk_index = spk_index + 1;
%     end 
    
% ---------------------------------------------------------------------------------

    if isempty(SpikeTrainH)
        SpikeTrainClusterSec = [];
    elseif clusterCount > length(SpikeTrainH) 
        SpikeTrainCluster = kmeans(SpikeTrainH', 1, 'Replicates',5, 'MaxIter', 500);
        SpikeTrainClusterSec = zeros(1, length(SpikeTrainCluster));
        SpikeTrainClusterSec(SpikeTrainCluster == 1) = clusterCount;
    else
        SpikeTrainCluster = kmeans(SpikeTrainH', clusterCount, 'Replicates',5, 'MaxIter', 500);
        
        clusterMaxValue = [];
        for i = 1:clusterCount
            clusterMaxValue(i) = max(SpikeTrainH(SpikeTrainCluster == i));
        end

        [~, cluster_sort_index] = sort(clusterMaxValue);
        SpikeTrainClusterSec = zeros(1, length(SpikeTrainCluster));
        for i = 1:clusterCount
            SpikeTrainClusterSec(SpikeTrainCluster == cluster_sort_index(i)) = i;
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
    
    for i = 1:clusterCount
        plot(max_index(SpikeTrainClusterSec == i), SpikeTrainH(SpikeTrainClusterSec == i), 'o');
    end
    xlim([1, size(traceSig, 1)]);
    ylim([-1, 5]);
    
    sb2 = subplot(8, 1, 8:8);
    imagesc(traceSig');
    colormap(jet);
    caxis(sb1.YLim);
       
    linkaxes([sb1, sb2], 'x');
    
    mysave(f, [outputpath, '\activity_averagePksHistByMLSpike_', num2str(activityIndex)]);
    
    comb_activity = fit;

    for i = 1:length(SpikeTrainPks)
        comb_activity(SpikeTrainStart(i):SpikeTrainEnd(i)) = traceSig(SpikeTrainStart(i):SpikeTrainEnd(i));
    end
end