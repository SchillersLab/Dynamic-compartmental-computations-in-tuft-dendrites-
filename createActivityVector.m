function createActivityVector(clusterType, curEventPearson, selectedROI, outputpathCurr, eventsTable, roiDistance, selectedROISplitDepth1, selectedROISplitDepth3)
    vectR = {};
    index = 1;
    Names = {};
    for i = 1:length(selectedROI)
        for k = (i+1):length(selectedROI)
            for e = 1:size(curEventPearson, 1)
                vectR(index, 1) = {sprintf('%s_%s', selectedROI{i}, selectedROI{k})};
                
                vectR(index, 2) = {curEventPearson(e, i, k)};
                vectR(index, 3) = eventsTable.event_name(e);
                
                vectR(index, 4) = {eventsTable.clusterByH(e)};
                vectR(index, 5) = {eventsTable.clusterByRoiPrecantage(e)};
                vectR(index, 6) = {roiDistance(i, k)};
                vectR(index, 7) = {selectedROISplitDepth1(i) == selectedROISplitDepth1(k) && selectedROISplitDepth1(i) ~= -1};
                vectR(index, 8) = {selectedROISplitDepth3(i) == selectedROISplitDepth3(k) && selectedROISplitDepth3(i) ~= -1};
                vectR(index, 9) = {sscanf(eventsTable.event_name{e}, 'event_%d')};
                vectR(index, 10) = {eventsTable.MantelTest(e)};
                
                for z = 1:(length(eventsTable.behave{e}))
                    lc = find(strcmp(Names, eventsTable.behave{e}{z}));
                    
                    if isempty(lc)
                        Names(end+1) = eventsTable.behave{e}(z);
                        vectR(index, end+1) = {1};
                    else
                        vectR(index, lc+10) = {1};
                    end                
                end
                
                index = index + 1;
            end
        end
    end
    
    tr = cell2table(vectR);
    tr.Properties.VariableNames= [{'Name', 'Pearson', 'EventName', 'ClusterByH', 'ClusterByP', 'Distance',...
        'Depth1', 'Depth2', 'eventNum', 'Pval_Mantel'}, Names];
    writetable(tr, fullfile(outputpathCurr, clusterType, 'matrixForGLM.csv'));
end