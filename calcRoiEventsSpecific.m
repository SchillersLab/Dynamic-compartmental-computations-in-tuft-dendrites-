function calcRoiEventsSpecific(tableEventsBehave, roiActivityNames, outputpath, behaveType)
    resultsCell(1,1:length(behaveType)) = behaveType;
    
    resultsCell(3, 1) = {'all clusters'};
    resultsCell(2, 1) = {'ROIs Name'};

    clustres = unique(tableEventsBehave.clusterByH);
    for i_c = 1:length(clustres)
        resultsCell(3+i_c, 1) = {sprintf('cluster_%d', clustres(i_c))};
    end
    
    eventsROIMat = reshape(cell2mat(tableEventsBehave.roisEvent), length(roiActivityNames),size(tableEventsBehave, 1));
    
    for roiIndex = 1:length(roiActivityNames)
        precentageOfEvents = sum(eventsROIMat(roiIndex, :)) ./ size(tableEventsBehave, 1);
        
        resultsCell(2, 1+roiIndex) = roiActivityNames(roiIndex);
        resultsCell(3, 1+roiIndex) = {precentageOfEvents};
        
        for i_c = 1:length(clustres)
            locationCluster = find(tableEventsBehave.clusterByH == clustres(i_c));
            precentageOfEventscluster = sum(eventsROIMat(roiIndex, locationCluster)) ./ size(locationCluster, 1);
            resultsCell(3+i_c,  1+roiIndex) = {precentageOfEventscluster};
        end
    end
    
    writetable(cell2table(resultsCell), fullfile(outputpath ,'ROIsEventsPrecentageForBehave.csv'));
end