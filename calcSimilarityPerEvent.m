function activtiyMatrixSum = calcSimilarityPerEvent(allEventsTable, i_cluster, selectedROI, roiActivity_comb, clusteringV)
    if i_cluster == 0
        currentEvents = allEventsTable;
    elseif i_cluster == -1
        currentEvents = allEventsTable(clusteringV ~= 1, :);
    else    
        currentEvents = allEventsTable(clusteringV == i_cluster, :);
    end
    
    activtiyMatrixSum = zeros(length(selectedROI), length(selectedROI));
    
    for event_i = 1:size(currentEvents, 1)
        activtiyMatrix = nan(length(selectedROI), length(selectedROI));
        locationE = currentEvents.start(event_i):currentEvents.pks(event_i);
        
        if length(locationE) == 1
            locationE = currentEvents.start(event_i):(currentEvents.pks(event_i) + 1);
        end
        
        for roi_i = 1:length(selectedROI)
            for roi_k = roi_i:length(selectedROI)
                corrEvents = corr([roiActivity_comb(locationE, roi_i),...
                roiActivity_comb(locationE, roi_k)], 'type', 'Pearson');
                
                activtiyMatrix(roi_i, roi_k) = corrEvents(1,2);
                activtiyMatrix(roi_k, roi_i) = corrEvents(1,2);
            end
        end
        locationE = [];
        activtiyMatrixSum = activtiyMatrixSum + activtiyMatrix;  
    end
    
    activtiyMatrixSum = activtiyMatrixSum / size(currentEvents, 1);
end