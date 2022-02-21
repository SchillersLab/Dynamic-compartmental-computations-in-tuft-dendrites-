function saveROIBranchingIndexAsLabel(eventWin, gRoi, selectedROISplitDepth1, selectedROI, outputpath, depthNumber, allEventsTable, roiActivity_comb, tr_frame_count)
    classesMDepth1 = unique(selectedROISplitDepth1);
    eventStr = [];
    selectedROISplitDepthToSave1 = selectedROISplitDepth1;
    for indexC = 1:length(classesMDepth1)
        if classesMDepth1(indexC) == -1
            eventStr = [eventStr 'ND' '_'];
            labelsLUT(indexC) = {'ND'};
        
        else
            eventStr = [eventStr gRoi.Nodes(classesMDepth1(indexC),:).Name{1} '_'];
            labelsLUT(indexC) = gRoi.Nodes(classesMDepth1(indexC),:).Name(1);
        end
        
        cls(indexC, :) = getTreeColor('within', indexC, true);
        selectedROISplitDepthToSave1(selectedROISplitDepth1 == classesMDepth1(indexC)) = indexC;
    end
    
    roiTableLabelDepth1.roiNames =  selectedROI;
    roiTableLabelDepth1.labelsLUT = labelsLUT;
    roiTableLabelDepth1.eventsStr = eventStr;    
    roiTableLabelDepth1.cls = cls;
    roiTableLabelDepth1.labels = selectedROISplitDepthToSave1;
    
    for events_index = 1:size(allEventsTable, 1)
        event_curr_loc = (allEventsTable.pks(events_index)-eventWin):(allEventsTable.pks(events_index)+eventWin);
        
        if sum(event_curr_loc < 0) > 0
            event_curr_loc = 1:(eventWin*2);
        end
        
        if sum(event_curr_loc > size(roiActivity_comb, 1)) > 0
            event_curr_loc = (size(roiActivity_comb, 1) - (eventWin*2)):size(roiActivity_comb, 1);
        end
        
        for roiIndex = 1:length(selectedROI)
            roiTableLabelDepth1.activity.dataEvents(roiIndex, 1:length(event_curr_loc), events_index) = roiActivity_comb(event_curr_loc,roiIndex);
        end
        
        roiTableLabelDepth1.activity.labels(events_index) = allEventsTable.clusterByH(events_index);
    end
    
    for trialIndex = 1:(size(roiActivity_comb, 1) / tr_frame_count)
         for roiIndex = 1:length(selectedROI)
            roiTableLabelDepth1.activity.dataTrials(roiIndex, 1:tr_frame_count, trialIndex) = roiActivity_comb(((trialIndex - 1) * tr_frame_count + 1):(trialIndex * tr_frame_count),roiIndex);
         end  
    end    
    
    save([outputpath, '\structuralTreeLabels_depth', num2str(depthNumber) '.mat'], 'roiTableLabelDepth1');
end

%     saveROIBranchingIndexAsLabel(eventWin, gRoi, selectedROISplitDepth1, selectedROI, outputpath, firstDepthCompare, allEventsTable, roiActivity_comb, tr_frame_count);
%     saveROIBranchingIndexAsLabel(eventWin, gRoi, selectedROISplitDepth3, selectedROI, outputpath, secDepthCompare, allEventsTable, roiActivity_comb, tr_frame_count);
