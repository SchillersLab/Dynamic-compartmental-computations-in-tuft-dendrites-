function [windowFULL, windowToPeak, loc, activityClusterValue, clusterCount] = calcActivityEventsWindowsAndPeaksType2(roiActivity, outputpath, threshold_std, PeakWidth, clusterPrecentage, histBinWidth, samplingRate, tr_frame_count, aV)
    frameNumForSec = 0.5 / (1 / samplingRate);

    DeconvType              = 1;       % type of deconvolution to use
    DeconvFiltTau           = .8;       % filter slope in 1/sec (N.A. in Fast Rise)
    DeconvFiltDur           = .5;       % smoothing filter duration in sec
    DeconvFiltRiseTime      = .3;       % rise time in dF/F for peak detection
    DeconvFiltRiseAmp       = .1;       % rise value Max-Min in dF/F for peak detection
    DeconvRespMinWidth      = .2;       % min cell response dF/F duration in sec
    DeconvRespMaxWidth      = 2.5;        % max cell response dF/F duration in sec
    
    all_locationFull_start = [];
    
    all_locationFull_end = [];
    
    all_locationFull_pks = [];
    
    all_locationFull_H = [];
    
    all_locationFull_W = [];
    
%     [coeff,score,latent,tsquared,explained,mu] = pca(roiActivity);
    
    for i = 1:size(roiActivity, 2)
        
%         for k = 1:tr_frame_count:size(roiActivity, 1)
%             ac_curr = roiActivity(k:(k + tr_frame_count - 1), i);

            ac_curr = roiActivity(:, i);

%             [SpikeTrainStart, SpikeTrainEnd, SpikeTrainPks, SpikeTrainH] = calcRoiEventDetectorByOri(ac_curr, DeconvFiltDur, DeconvFiltRiseTime,DeconvRespMinWidth, DeconvRespMaxWidth, 1 / samplingRate);
%             [SpikeTrainStart, SpikeTrainEnd, SpikeTrainPks, SpikeTrainH, SpikeTrainW] = calcRoiEventDetectorBySlop(ac_curr, DeconvRespMinWidth, DeconvRespMaxWidth, 1 / samplingRate, DeconvFiltDur);
            [SpikeTrainStart, SpikeTrainEnd, SpikeTrainPks, SpikeTrainH] = calcRoiEventDetectorByMLSpike(ac_curr, 1 / samplingRate, tr_frame_count, aV(i));

%             SpikeTrainStart = SpikeTrainStart + k -1;
%             SpikeTrainEnd = SpikeTrainEnd + k -1;
%             SpikeTrainPks = SpikeTrainPks + k -1;
%             
    
    %             
            all_locationFull_start(end+1:end+length(SpikeTrainStart)) = SpikeTrainStart;
            all_locationFull_end(end+1:end+length(SpikeTrainStart)) = SpikeTrainEnd;
            all_locationFull_pks(end+1:end+length(SpikeTrainStart)) = SpikeTrainPks;
            all_locationFull_H(end+1:end+length(SpikeTrainStart)) = SpikeTrainH;
%             all_locationFull_W(end+1:end+length(SpikeTrainStart)) = SpikeTrainW;

            %         end
       
    end

    [all_locationFull_pksSort, indexSortLoc] = sort(all_locationFull_pks);
    all_locationFull_endSort = all_locationFull_end(indexSortLoc);
    all_locationFull_startSort = all_locationFull_start(indexSortLoc);
    all_locationFull_HSort = all_locationFull_H(indexSortLoc);
%     all_locationFull_WSort = all_locationFull_W(indexSortLoc);
    
    for index = 1:length(all_locationFull_pksSort) - 1
        if abs(all_locationFull_pksSort(index+1) - all_locationFull_pksSort(index)) <= samplingRate/4
            
            all_locationFull_pksSort(index + 1) = round((all_locationFull_pksSort(index + 1) + all_locationFull_pksSort(index)) / 2);
            all_locationFull_startSort(index + 1) = round((all_locationFull_startSort(index + 1) + all_locationFull_startSort(index)) / 2);
            all_locationFull_endSort(index + 1) = round((all_locationFull_endSort(index + 1) + all_locationFull_endSort(index)) / 2);
            all_locationFull_HSort(index + 1) = (all_locationFull_HSort(index + 1) + all_locationFull_HSort(index)) / 2;
%             all_locationFull_WSort(index + 1) = (all_locationFull_WSort(index + 1) + all_locationFull_WSort(index)) / 2;
            
            all_locationFull_pksSort(index) = -1;
            all_locationFull_startSort(index) = -1;
            all_locationFull_endSort(index) = -1;
            all_locationFull_HSort(index) = -1;
%             all_locationFull_WSort(index) = -1;
%         else
%             if index ~= 1 && all_locationFull_pksSort(index-1) ~= -1
%                 all_locationFull_pksSort(index) = -1;
%                 all_locationFull_startSort(index) = -1;
%                 all_locationFull_endSort(index) = -1;
%                 all_locationFull_HSort(index) = -1;
%             end
        end
    end
    
    all_locationFull_pksSort(all_locationFull_pksSort == -1) = [];
    all_locationFull_startSort(all_locationFull_startSort == -1) = [];
    all_locationFull_endSort(all_locationFull_endSort == -1) = [];
    all_locationFull_HSort(all_locationFull_HSort == -1) = [];
%     all_locationFull_WSort(all_locationFull_WSort == -1) = [];
  
    activity_average = mean(roiActivity,2);     
    loc = all_locationFull_pksSort;
    loc2 = [];
    pks = [];
    for index = 1:length(all_locationFull_startSort)
        [maxV, maxL] = max(activity_average(all_locationFull_startSort(index):all_locationFull_endSort(index)));
        pks(index) = maxV;
        loc2(index) = maxL + all_locationFull_startSort(index) - 1;
    end
    
    fig = figure;
    hold on;
    T2 = histogram(pks', 'BinWidth', histBinWidth, 'Normalization', 'pdf'); 
    pd2 = fitdist(pks','Kernel');
    y2 = pdf(pd2,T2.BinLimits(1):0.01:T2.BinLimits(2));
    plot(T2.BinLimits(1):0.01:T2.BinLimits(2), y2);
    
    mysave(fig, [outputpath, '\activity_averagePksHist']);
         
     for i = 1:15
        [idx(i, :),~,sumdT] = kmeans(pks', i, 'Replicates',5, 'MaxIter', 500);   
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
   
    for i = 1:clusterCount
        plot(loc2(activityClusterValue == i), pks(activityClusterValue == i),'o', 'DisplayName', ['cluster', num2str(i)]);
        fprintf(fileID,'%d Cluster - %d events max value %d, min value %d \n', sum(activityClusterValue == i), i,  max(pks(activityClusterValue == i)), min(pks(activityClusterValue == i)));
    end
    
    legend('show');
    mysave(fig, [outputpath, '\SelectedEventsForROIS']);
    fclose(fileID);
    
    windowFULL = [all_locationFull_startSort, all_locationFull_endSort];    
    windowToPeak = [all_locationFull_startSort, all_locationFull_pksSort];    
    eventsWindowsActivity_fromPeak = [all_locationFull_pksSort, all_locationFull_endSort];          
end 