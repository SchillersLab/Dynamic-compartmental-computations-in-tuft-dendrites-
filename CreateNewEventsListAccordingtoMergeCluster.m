function CreateNewEventsListAccordingtoMergeCluster
    outputpath = "";
    roiActivityCombList = {'', ''};
    clusterCount = 4;
    
    [allEventsTableComb, filesActivity.nume1.roiActivity_comb] = load([outputpath, '\roiActivity_comb.mat'], 'allEventsTable', 'roiActivity_comb');
    allEventsTableComb.fileNumber = ones(size(allEventsTableComb, 1), 1);
    
    for i = 2:length(roiActivityCombList)
        [allEventsTable, filesActivity.(sprintf('nume%d', i)).roiActivity_comb] = load([outputpath, '\roiActivity_comb.mat'], 'allEventsTable', 'roiActivity_comb');
        allEventsTable.fileNumber = ones(size(allEventsTable, 1), 1) * i;
        
        allEventsTableComb(end+1:end+size(allEventsTable, 1), :) = allEventsTable;
    end
    
    SpikeTrainClusterSecByH = getClusterForActivity(allEventsTableComb.H, clusterCount);
    SpikeTrainClusterSecByP = getClusterForActivity(allEventsTableComb.roiPrecantage, clusterCount);
    
    allEventsTableComb.clusterByH = SpikeTrainClusterSecByH;
    allEventsTableComb.clusterByRoiPrecantage = SpikeTrainClusterSecByP;
    
    for i = 1:length(roiActivityCombList)
         allEventsTable = allEventsTableComb(allEventsTableComb.fileNumber == i,:);
         roiActivity_comb = filesActivity.(sprintf('nume%d', i)).roiActivity_comb;
         save([outputpath, '\roiActivity_comb', i,'.mat'], 'allEventsTable', 'roiActivity_comb');
    end
end