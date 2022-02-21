function [windowFULL, windowToPeak, loc, activityClusterValue, clusterCount] = calcActivityEventsWindowsAndPeaks(roiActivity, outputpath, threshold_std, PeakWidth, clusterPrecentage, histBinWidth, samplingRate)
    frameNumForSec = 0.5 / (1 / samplingRate);
% 
%     DeconvType              = 1;       % type of deconvolution to use
%     DeconvFiltTau           = .8;       % filter slope in 1/sec (N.A. in Fast Rise)
%     DeconvFiltDur           = .2;       % smoothing filter duration in sec
%     DeconvFiltRiseTime      = .1;       % rise time in dF/F for peak detection
%     DeconvFiltRiseAmp       = .2;       % rise value Max-Min in dF/F for peak detection
%     DeconvRespMinWidth      = .2;       % min cell response dF/F duration in sec
%     DeconvRespMaxWidth      = 2;        % max cell response dF/F duration in sec
%     
%     all_locationFull_start = [];
%     
%     all_locationFull_end = [];
%     
%     all_locationFull_pks = [];
%     
%     for i = 1:size(roiActivity, 2)
%         calcRoiEventDetector(roiActivity(:, i), DeconvType, DeconvFiltTau, DeconvFiltDur, DeconvFiltRiseTime, DeconvFiltRiseAmp, DeconvRespMinWidth, DeconvRespMaxWidth, samplingRate)
%         
%         windowFULL = zeros(length(loc_cur), 2);    
%         windowToPeak = zeros(length(loc_cur), 2);    
%        
%         for index = 1:length(loc_cur)
%             windowFULL = [loc_cur - round(w_cur / 2), loc_cur + round(w_cur / 2)];
%             windowToPeak = [loc_cur - round(w_cur / 2), loc_cur];
%             
%         end
%         
%         deletedIndexList = [];
%         for index = 1:size(windowToPeak, 1)
%             if (windowToPeak(index, 2) - windowToPeak(index,1) < frameNumForSec)
%                 deletedIndexList(end+1) = index;
%             end
%         end
%         
%         windowFULL(deletedIndexList, :) = [];
%         windowToPeak(deletedIndexList, :) = [];
%         
%         all_locationFull_start(end+1 : end + size(windowToPeak, 1)) = windowToPeak(:, 1);
%         all_locationFull_end(end+1 : end + size(windowToPeak, 1)) = windowFULL(:, 2);
%         all_locationFull_pks(end+1 : end + size(windowToPeak, 1)) = windowToPeak(:, 2);
%     end
% 
%     [~, indexSortLoc] = sort(all_locationFull_start);
%     
%     combLoc = indexSortLoc(1);
%     
%     for indexL = indexSortLoc(2:end)
%         if all_locationFull_start(combLoc(end)) ~= all_locationFull_start(indexL) &&...
%            (all_locationFull_start(indexL) - all_locationFull_start(combLoc(end))) > frameNumForSec / 2
%             combLoc(end + 1) = indexL;
%         end
%     end
%     
%     
    activity_average = mean(roiActivity,2);     
%     XDEN = activity_average;
    XDEN = wdenoise(activity_average);   
 
    fig = figure;
    hold on;
    
    T = histogram(XDEN, 'BinWidth', histBinWidth, 'Normalization', 'pdf');
    pd = fitdist(XDEN,'Kernel');
    y = pdf(pd,T.BinLimits(1):0.01:T.BinLimits(2));
    plot(T.BinLimits(1):0.01:T.BinLimits(2), y);
    
    maxBinlocation = find(T.Values == max(T.Values), 1);
    stdA = std(XDEN(XDEN <= T.BinEdges(maxBinlocation)));
    threshold = T.BinEdges(maxBinlocation) + threshold_std * stdA;
    baseLine = T.BinEdges(maxBinlocation);
    mysave(fig, [outputpath, '\activity_averageHist']);
 
    
%     baseLine = iqr(activity_average);
%     threshold = threshold_std*baseLine;
    [pks,loc,w_loc,~]  = findpeaks(XDEN,'MinPeakHeight',threshold, 'MinPeakWidth', PeakWidth, 'MinPeakDistance', PeakWidth);
      
    fig = figure;
    hold on;
    T2 = histogram(pks, 'BinWidth', histBinWidth, 'Normalization', 'pdf'); 
    pd2 = fitdist(pks,'Kernel');
    y2 = pdf(pd2,T2.BinLimits(1):0.01:T2.BinLimits(2));
    plot(T2.BinLimits(1):0.01:T2.BinLimits(2), y2);
    
    mysave(fig, [outputpath, '\activity_averagePksHist']);
         
     for i = 1:15
        [idx(i, :),~,sumdT] = kmeans(pks, i, 'Replicates',5, 'MaxIter', 500);   
        sumd(i) = sum(sumdT);
     end
     
     if isempty(sumd > 1)
        sumd = sumd .* 100;
     end
     
    clusterCount = find(sumd < clusterPrecentage, 1);

    activityClusterValue = idx(clusterCount, :);
     
    fig = figure;
    hold on;
    plot(sumd);
    plot(clusterCount, sumd(clusterCount), 'ro');  
    mysave(fig, [outputpath, '\clusteringSelection']);   
    
    fileID = fopen(fullfile(outputpath, 'EventClusterResults.txt'),'w');
    
    fig = figure;
    hold on;
    plot(1:length(activity_average), activity_average, 'DisplayName', 'activity');
    plot(1:length(XDEN), XDEN, 'DisplayName', 'denois');

    for i = 1:clusterCount
        plot(loc(activityClusterValue == i), pks(activityClusterValue == i),'o','MarkerSize',12, 'DisplayName', ['cluster', num2str(i)]);
        fprintf(fileID,'%d Cluster - %d events max value %d, min value %d \n', sum(activityClusterValue == i), i,  max(pks(activityClusterValue == i)), min(pks(activityClusterValue == i)));
    end
    
    plot(1:length(activity_average), threshold * ones(length(activity_average),1), 'DisplayName', 'threshold');
    legend('show');
    mysave(fig, [outputpath, '\SelectedEventsForROIS']);
    fclose(fileID);
    
     
    fig = figure;
    hold on;
    plot(1:length(XDEN), XDEN, 'DisplayName', 'denois');
    
    for i = 1:clusterCount
        plot(loc(activityClusterValue == i), pks(activityClusterValue == i),'o','MarkerSize',12, 'DisplayName', ['cluster', num2str(i)]);
     end
    
    plot(1:length(activity_average), threshold * ones(length(activity_average),1), 'DisplayName', 'threshold');
    legend('show');
    mysave(fig, [outputpath, '\SelectedEventsForROISOnlyDenois']);
    
    windowFULL = zeros(length(loc), 2);    
    windowToPeak = zeros(length(loc), 2);    
    eventsWindowsActivity_fromPeak = zeros(length(loc), 2);    
%     
%     for index = 1:length(loc)
%         belowBaseLine_b = find(XDEN(loc(index):-1:1) <= baseLine,1);
%         if isempty(belowBaseLine_b)
%             [~, belowBaseLine_b] = min(XDEN(loc(index):-1:1));
%         end
%         
%         windowFULL(index, 1) = loc(index) - belowBaseLine_b + 1;
%         
%         belowBaseLine = find(XDEN(loc(index):1:length(XDEN)) <= baseLine,1);
%         if isempty(belowBaseLine)
%             [~, belowBaseLine] = min(XDEN(loc(index):1:length(XDEN)));
%         end
%         
%         windowFULL(index, 2) = loc(index) + belowBaseLine - 1;
%         
%         if (index > 1 && (windowFULL(index, 1) <= loc(index - 1, 1)))
%             act = XDEN(loc(index-1):1:loc(index));
%             TF = loc(index-1) + find(act == min(act)) - 1;
%             windowFULL(index, 1) = TF;
%             windowFULL(index - 1, 2) = TF;  
%             eventsWindowsActivity_fromPeak(index - 1, 2) = TF;
%         end
%         
%         windowToPeak(index, 1) = windowFULL(index, 1);
%         windowToPeak(index, 2) = loc(index);
%         
%         eventsWindowsActivity_fromPeak(index, 1) = loc(index);
%         eventsWindowsActivity_fromPeak(index, 2) = windowFULL(index, 2); 
%     end  
end 