function [SpikeTrainStart, SpikeTrainEnd, SpikeTrainPks, SpikeTrainH, SpikeTrainW] = calcRoiEventDetectorBySlop(dataCurROI, DeconvRespMinWidth, DeconvRespMaxWidth, ImageSampleTime, DeconvFiltDur)       
    %%%%%
    % Define Filter
    %%%%%
    sampFreq        = 1/ImageSampleTime;
    
    %MinSpikeWidth       = ceil(filtTau);     % samples for average spike
    respMinWidth        = ceil(DeconvRespMinWidth * sampFreq); % minimal width
    respMaxWidth        = ceil(DeconvRespMaxWidth * sampFreq); % max width

    % Processing works on columns
    traceSig            = dataCurROI;
    nR = size(dataCurROI, 1);
    % smooth
    
    filtDur         = DeconvFiltDur * sampFreq;      % filter duration in sec
    %filtTau     = DeconvFiltTau * sampFreq;     % filter slope in 1/sec
    filtLenH        = ceil(filtDur/2);
    filtLen         = filtLenH*2;

    % preobje for deconvolution - freq domain
    %filtSmooth  = ones(filtLen,1);
    filtSmooth          = hamming(filtLen);
    %filtSmooth  = gausswin(filtLen);
    filtSmooth          = filtSmooth./sum(filtSmooth);
    
    dFFFilt2             = filtfilt(filtSmooth,1,traceSig);
    
    
    dFFFilt = wdenoise(traceSig);
    noiseSig = traceSig - dFFFilt;
    baseLine = iqr(noiseSig);
    stdNoise = std(noiseSig);
    
    MaxMinThr = baseLine + stdNoise * 1;
    
    t = 1;
    v_calc = (dFFFilt2(2:t:(end)) - dFFFilt2(1:t:(end-1))) ./ t;
    v_calc = repelem(v_calc,t,1);
    v_calc = [v_calc(1) ;v_calc];
    v_slop = v_calc > 0;
    
    res = [0; find(v_slop(1:end-1) ~= v_slop(2:end)); nR];
    
    for k = 1:length(res)-1
       if v_slop(res(k) + 1) == 0 & (res(k + 1) - res(k) <= 2)
            slop_value = 1;
            if (k-1>0)
                slop_value = v_slop(res(k-1) + 1) ;    
            end
            
            v_slop((res(k) + 1):res(k+1)) = slop_value;
        end
    end
    
   res = [0; find(v_slop(1:end-1) ~= v_slop(2:end)); nR];
   
    
    for k = 1:length(res)-1
        if v_slop(res(k) + 1) == 1 & (res(k + 1) - res(k) <= 2)
            slop_value = 0;
            
            if (k-1>0)
                slop_value = v_slop(res(k-1) + 1) ;    
            end
            
            v_slop((res(k) + 1):res(k+1)) = slop_value;
        end
    end
    
    
    res = [0; find(v_slop(1:end-1) ~= v_slop(2:end)); nR];
    
    %     % create spike train
    SpikeTrainStart  = [];
    SpikeTrainEnd  = [];
    SpikeTrainPks  = [];
    SpikeTrainH  = [];
    SpikeTrainAngle  = [];
    SpikeTrainW  = [];
    
    dMaxMinInd = [];
    
    resIndex = 1;
    while (resIndex + 2 <= length(res))
        if (v_slop(res(resIndex) + 1) == 1) 
            resIndexSec = resIndex + 2;
                
            isAboveThreshold = find(dFFFilt2((res(resIndex) + 1):res(resIndexSec)) > MaxMinThr);       
            
            if ~isempty(isAboveThreshold) && length(isAboveThreshold) >= 1
                [SpikeH, SpikeLoc] = max(dFFFilt2((res(resIndex) + 1):res(resIndex + 1)));
                SpikeLoc = SpikeLoc + res(resIndex);
                
                end_pos = find((dFFFilt2(SpikeLoc:res(resIndexSec))) <= baseLine, 1);
                if isempty(end_pos)
                    end_pos = res(resIndexSec);
                else
                    end_pos = end_pos + SpikeLoc - 1;
                end
                
                start_pos = find((dFFFilt2((res(resIndex) + 1):SpikeLoc)) <= baseLine, 1, 'last');
                if isempty(start_pos)
                    start_pos = (res(resIndex) + 1);
                else
                    start_pos = start_pos + res(resIndex);
                end
                
                
                widthSpike = find(dFFFilt2(start_pos:end_pos) >= SpikeH * 0.5);
                widthSpike = length(widthSpike);  

                if (SpikeLoc - start_pos == 0)
                    resIndex = resIndex + 2;
                    continue;
                end
                
                slop_rise_value = (dFFFilt2(SpikeLoc) - dFFFilt2(start_pos)) ./ (SpikeLoc - start_pos);
%                 slop_rise_angle = (180 / pi ) * atan(slop_rise_value);

                slop_down_value = abs((dFFFilt2(end_pos) - dFFFilt2(SpikeLoc)) ./ (end_pos - SpikeLoc));
%                 slop_down_angle = (180 / pi ) * atan(slop_down_value);

                if ~isnan(slop_rise_value) && ~isnan(slop_down_value) && abs(slop_down_value) <= abs(slop_rise_value) && widthSpike >= 3 
                    currWidthToPeak = SpikeLoc - start_pos;
%                     currWidthFPeak = end_pos - SpikeLoc;
                    if (currWidthToPeak >= respMinWidth && currWidthToPeak < respMaxWidth) 
%                             (currWidthFPeak >= respMinWidth && currWidthFPeak < respMaxWidth) 
                        SpikeTrainStart(end + 1) = start_pos;
                        SpikeTrainEnd(end + 1) = end_pos;
                        SpikeTrainPks(end + 1) = SpikeLoc;
                        SpikeTrainH(end + 1) = SpikeH;
%                         SpikeTrainAngle(end + 1, :) = [slop_rise_angle, slop_down_angle];
                        SpikeTrainW(end + 1) = widthSpike;
                    end
                end
            end
            
            resIndex = resIndex + 2;
        else            
            resIndex = resIndex + 1;
        end
    end
    
%     for index = 1:length(SpikeTrainPks) - 1
%         if SpikeTrainPks(index+1) - SpikeTrainPks(index) <= sampFreq/2
%             if (dFFFilt2(SpikeTrainPks(index+1)) < dFFFilt2(SpikeTrainPks(index)))
%                 SpikeTrainPks(index+1) = SpikeTrainPks(index);
%                 SpikeTrainStart(index + 1) = SpikeTrainStart(index);
%                 SpikeTrainH(index + 1) = SpikeTrainH(index);
%             else
%                 SpikeTrainStart(index + 1) = SpikeTrainStart(index);
%             end
%             
%             SpikeTrainPks(index) = -1;
%             SpikeTrainEnd(index) = -1;
%             SpikeTrainStart(index) = -1;
%             SpikeTrainH(index) = -1;
% %             SpikeTrainAngle(index, :) = [-1, -1];
%             SpikeTrainW(index) = -1;
%         end
%     end
%     
%     SpikeTrainStart(SpikeTrainPks == -1) = [];
%     SpikeTrainEnd(SpikeTrainPks == -1) = [];    
%     SpikeTrainH(SpikeTrainPks == -1) = [];
%     SpikeTrainW(SpikeTrainPks == -1) = [];
% %     SpikeTrainAngle(SpikeTrainPks == -1) = [];
%     SpikeTrainPks(SpikeTrainPks == -1) = [];
%     
    figure;
    hold on;
    plot(traceSig)
    plot(dFFFilt2)
    plot(SpikeTrainPks, dFFFilt2(SpikeTrainPks), '*r');
    plot(SpikeTrainEnd, dFFFilt2(SpikeTrainEnd), '*b');
    plot(SpikeTrainStart, dFFFilt2(SpikeTrainStart), '*g');
    plot(1:nR,  MaxMinThr * ones(nR,1), 'DisplayName', 'threshold'); 
end