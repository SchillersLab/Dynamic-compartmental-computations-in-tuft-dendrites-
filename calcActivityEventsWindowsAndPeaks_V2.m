function [all_locationFull_start, all_locationFull_end, all_locationFull_pks, all_locationFull_cluster, roiActivity_comb] = calcActivityEventsWindowsAndPeaks_V2(roiActivity, outputpath, clusterCount, samplingRate, tr_frame_count, aV, roiActivityNames, sigmaChangeValue)
    
    all_locationFull_start = [];
    
    all_locationFull_end = [];
    
    all_locationFull_pks = [];
    
    all_locationFull_H = [];
    
    all_locationFull_cluster = [];
    
    for i = 1:size(roiActivity, 2)
        ac_curr = roiActivity(:, i);
        [all_locationFull_start{i}, all_locationFull_end{i}, all_locationFull_pks{i}, all_locationFull_H{i}, all_locationFull_cluster{i}, roiActivity_comb(:, i)] = calcRoiEventDetectorByMLSpike(ac_curr, 1 / samplingRate, tr_frame_count, aV(i), outputpath, i, clusterCount, roiActivityNames(i), sigmaChangeValue(i));
    end

%     idx = [];
%     sumd = [];
%     for k = 1:15
%         [idx(k, :),~,sumdT] = kmeans(all_locationFull_H_forCluster', k, 'Replicates',5, 'MaxIter', 500);   
%         sumd(k) = sum(sumdT);
%     end
%     
%     sumd = (sumd ./ sum(sumd)) .* 100;
% 
% %     if isempty(sumd > 1)
% %     sumd = sumd .* 100;
% %     end
% 
%     clusterCount = find(sumd < clusterPrecentage, 1);
% 
%     fig = figure;
%     hold on;
%     plot(sumd);
%     plot(clusterCount, sumd(clusterCount), 'ro');  
%     mysave(fig, [outputpath, '\clusteringSelection']);    
%     
%     fileID = fopen(fullfile(outputpath, 'EventClusterResults.txt'),'w');
%     
%     fig = figure;
%     hold on;
%    
%     for i = 1:clusterCount
%         plot(all_locationFull_H_forCluster(idx(clusterCount, :) == i), '*', 'DisplayName', ['cluster', num2str(i)]);
%         fprintf(fileID,'%d Cluster - %d events max value %d, min value %d \n', sum(idx(clusterCount, :) == i), i,  max(all_locationFull_H_forCluster(idx(clusterCount, :) == i)), min(all_locationFull_H_forCluster(idx(clusterCount, :) == i)));
%     end
%     
%     legend('show');
%     mysave(fig, [outputpath, '\SelectedEventsForROIS']);
%     fclose(fileID);
%     
%     startPoint = 1;
%     for i = 1:length(all_locationFull_H)
%         endPoint = startPoint + length(all_locationFull_H{i}) - 1; 
%         activityClusterValue{i} = idx(clusterCount, startPoint:endPoint);
%         
%         startPoint = startPoint + length(all_locationFull_H{i});
%     end
end