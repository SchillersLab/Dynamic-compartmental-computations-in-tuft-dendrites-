tableResults = tableResults2;

allR2Adj.c0 = [];
allR2Adj.c5 = [];
allR2Adj.c1 = [];
allR2Adj.c2 = [];
allR2Adj.c3 = [];
allR2Adj.c4 = [];
for i = 1:size(tableResults, 1)
    load([tableResults.RunLocation{i}, '\roiActivityRawData.mat']);
    load([tableResults.RunLocation{i}, '\..\..\..\EventsDetection\roiActivity_comb.mat'], 'allEventsTable', 'roiActivity_comb');
    
    classesRoi = unique(selectedROISplitDepth1);
    numberOfminRoiPerSideR = sum(selectedROISplitDepth1==classesRoi(1));
    numberOfminRoiPerSideL = sum(selectedROISplitDepth1==classesRoi(2));
    
    if numberOfminRoiPerSideL < 8 | numberOfminRoiPerSideR < 8
        continue;
    end
    
    eventsIndexLocationCluster = [];
    eventsIndexLocation = [];
    for j = 1:size(allEventsTable, 1)
        eventsIndexLocation(end+1:end+(allEventsTable.pks(j)-allEventsTable.start(j))+1) = allEventsTable.start(j):allEventsTable.pks(j);
        eventsIndexLocationCluster(end+1:end+(allEventsTable.pks(j)-allEventsTable.start(j))+1) = ones(1, (allEventsTable.pks(j)-allEventsTable.start(j))+1) * allEventsTable.clusterByH(j);
    end

    for ci = 0:5        
        if ci == 0
            eventsIndexLocationCurrent = eventsIndexLocation;
        elseif ci == 5
            eventsIndexLocationCurrent = eventsIndexLocation(eventsIndexLocationCluster ~= 4);
        else
            eventsIndexLocationCurrent = eventsIndexLocation(eventsIndexLocationCluster == ci);
        end 
        
        currentVar = [];
        for j = 1:length(selectedROISplitDepth1)
            otherBranch = find(selectedROISplitDepth1 ~= selectedROISplitDepth1(j));
            mdlTest = fitglm(roiActivity_comb(eventsIndexLocationCurrent, otherBranch), roiActivity_comb(eventsIndexLocationCurrent, j));
            currentVar(end+1) = mdlTest.Rsquared.Adjusted;
        end

%         allR2Adj.(sprintf('c%d',ci))(end+1:end+length(currentVar)) = currentVar;
        allR2Adj.(sprintf('c%d',ci))(end+1) = mean(currentVar);
    end
end

test = 1;
