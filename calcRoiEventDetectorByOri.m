function [SpikeTrainStart, SpikeTrainEnd, SpikeTrainPks, SpikeTrainH] = calcRoiEventDetectorByOri(dataCurROI, DeconvFiltDur, DeconvFiltRiseTime, DeconvRespMinWidth, DeconvRespMaxWidth, ImageSampleTime)       
    %%%%%
    % Define Filter
    %%%%%
    sampFreq        = 1/ImageSampleTime;
    filtDur         = DeconvFiltDur * sampFreq;      % filter duration in sec
    %filtTau     = DeconvFiltTau * sampFreq;     % filter slope in 1/sec
    filtLenH        = ceil(filtDur/2);
    filtLen         = filtLenH*2;

    % preobje for deconvolution - freq domain
    %filtSmooth  = ones(filtLen,1);
    filtSmooth          = hamming(filtLen);
    %filtSmooth  = gausswin(filtLen);
    filtSmooth          = filtSmooth./sum(filtSmooth);

    % additional objams
%     MaxMinThr           = DeconvFiltRiseAmp;  % determine the difference in rize time
    
    
    %MinSpikeWidth       = ceil(filtTau);     % samples for average spike
    ArtifactTime        = 2;    % samples in the initial time to remove spikes
    SupWinSize          = ceil(DeconvFiltRiseTime * sampFreq);
    respMinWidth        = ceil(DeconvRespMinWidth * sampFreq); % minimal width
    respMaxWidth        = ceil(DeconvRespMaxWidth * sampFreq); % max width

    % Processing works on columns
    traceSig            = dataCurROI;
    nR = size(dataCurROI, 1);
    % smooth
    
    dFFFilt             = filtfilt(filtSmooth,1,traceSig);
    
    dFFFilt_dnoise = wdenoise(traceSig);
    noiseSig = traceSig - dFFFilt_dnoise;
    baseLine = iqr(noiseSig);
    stdNoise = std(noiseSig);
    
    MaxMinThr = baseLine + stdNoise * 2;
    

    %     % create spike train
    SpikeTrainStart  = [];
    SpikeTrainEnd  = [];
    SpikeTrainPks  = [];
    SpikeTrainH  = [];

%     -------------------------------------------------------------------
    dFFMinFix           = dFFFilt;
    dFFMaxFix           = dFFFilt*0;%[dFFFilt(SupWinSize:end); zeros(SupWinSize-1,1)+0 ];
    for n = 1:nR-SupWinSize
        dFFMaxFix(n) = max(dFFFilt(n:n+SupWinSize-1));
    end

    % find places with bif difference betwen min and max
    dMaxMin             = dFFMaxFix - dFFMinFix;            
    dMaxMin(1:ArtifactTime) = 0;  % if any artifact  
    dMaxMinInd          = find(dMaxMin > MaxMinThr & [dMaxMin(2:end); 0] < dMaxMin & dMaxMin >= [10; dMaxMin(1:end-1)]);

    s_num       = numel(dMaxMinInd);
    SpikeArea   = zeros(s_num,1); % contains response width
    SpikeWidth  = zeros(s_num,1); % contains response width
    SpikeHeight = zeros(s_num,1); % contains response width
    SpikeLoc = zeros(s_num,1);
    for s = 1:s_num
        first_point  = dMaxMinInd(s);
        
        last_point   = min(nR,first_point + 1); %min(nR,first_point + MinSpikeWidth);
        if s == s_num
            lastPos      = nR;
        else
            lastPos      = min(nR,dMaxMinInd(s+1) - 3);  % go to the next rise. small gap
        end
        
        localThr    = dFFFilt(first_point)*1.00; % max(dFFFilt(lastPos),dFFFilt(first_point))*1.1;
        
        while dFFFilt(last_point) >= localThr && last_point < lastPos
          last_point = last_point + 1;
        end
        
        % record only above minimal width
        currWidth       = (last_point - first_point);
        if currWidth >= respMinWidth && currWidth < respMaxWidth
            [SpikeHeight(s), SpikeLoc(s)] = max(dFFFilt(first_point:last_point));
%             SpikeHeight(s) = SpikeHeight(s) - dFFFilt(first_point);
            SpikeLoc(s) = SpikeLoc(s) + first_point - 1;
            SpikeArea(s)   = sum(dFFFilt(first_point:last_point) - dFFFilt(first_point)) ; % resolves when signal below 0
            SpikeWidth(s)  = currWidth;
            
            SpikeTrainStart(end + 1) = first_point;
            SpikeTrainEnd(end + 1) = last_point;
            SpikeTrainPks(end + 1) = SpikeLoc(s);
            SpikeTrainH(end + 1) = SpikeHeight(s);
         end
    end
    
    figure;
    hold on;
    plot(traceSig)
    plot(dFFFilt)
    plot(SpikeTrainPks, dFFFilt(SpikeTrainPks), 'o');
    plot(1:nR,  MaxMinThr * ones(nR,1), 'DisplayName', 'threshold'); 
end