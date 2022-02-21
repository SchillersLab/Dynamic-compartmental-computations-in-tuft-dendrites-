function plotEventsSperation(allEventsTable, roiActivity_comb, selectedROI, selectedROISplitDepth1, outputpath, gRoi, roiSortedByCluster, roiTreeDistanceMatrix, index_apical)
    classes = unique(selectedROISplitDepth1);
    classes(classes == -1) = [];
    eventSplit = zeros(size(allEventsTable, 1), length(classes));
    
    outputpath = [outputpath, '\EventsPlots\'];
    
    for ca_i = 1:size(allEventsTable, 1)        
        for cl = 1:length(classes)
            eventSplit(ca_i, cl) = sum(allEventsTable.roisEvent{ca_i}(selectedROISplitDepth1 == classes(cl))) / sum(selectedROISplitDepth1 == classes(cl));
        end
    end
    
    for cl = 1:(length(classes)-1)
        fig = figure;
        hold on;
        
        scatter(eventSplit(:, cl), eventSplit(:, cl+1), 'filled', 'SizeData', 20, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
        xlabel(sprintf('Precentage in %s', gRoi.Nodes.Name{classes(cl)}));
        ylabel(sprintf('Precentage in %s', gRoi.Nodes.Name{classes(cl+1)}));
        xlim([0,1]);
        ylim([0,1]);
        plot([0.5, 0.5], [0,1], '--k', 'LineWidth', 1.5);
        plot([0, 1], [0.5,0.5], '--k', 'LineWidth', 1.5);
        title('Precentage of Roi contributed to Ca Events');
        mysave(fig, [outputpath, 'EventsPlotByP']);
        
        snapnow;
        close all;
        
        
        eventsCl1Location = eventSplit(:, cl) > eventSplit(:, cl+1) & eventSplit(:, cl) > 0.2;
        eventsCl2Location = eventSplit(:, cl) < eventSplit(:, cl+1) & eventSplit(:, cl+1) > 0.2;
        eventsCl3Location = eventSplit(:, cl) > 0.5 & eventSplit(:, cl+1) > 0.5;
        eventsCl4Location = eventSplit(:, cl) > 0 & eventSplit(:, cl+1) <= 0;
        eventsCl5Location = eventSplit(:, cl) <= 0 & eventSplit(:, cl+1) > 0;
        
        sprintf('Ca events, subtree %s is more active', gRoi.Nodes.Name{classes(cl)})
        plotHMEvents(allEventsTable(eventsCl1Location, :), roiActivity_comb, selectedROI, [outputpath, '\',  gRoi.Nodes.Name{classes(cl)}], roiSortedByCluster, selectedROISplitDepth1, gRoi, roiTreeDistanceMatrix, index_apical);
        
        snapnow;
        close all;
        
        sprintf('Ca events, subtree %s is more active', gRoi.Nodes.Name{classes(cl+1)})
        plotHMEvents(allEventsTable(eventsCl2Location, :), roiActivity_comb, selectedROI, [outputpath, '\',  gRoi.Nodes.Name{classes(cl+1)}], roiSortedByCluster, selectedROISplitDepth1, gRoi, roiTreeDistanceMatrix, index_apical);
        
        snapnow;
        close all;
        
        sprintf('Ca events, subtree both roi"s above 50 precentage')
        plotHMEvents(allEventsTable(eventsCl3Location, :), roiActivity_comb, selectedROI, [outputpath, '\above50'], roiSortedByCluster, selectedROISplitDepth1, gRoi, roiTreeDistanceMatrix, index_apical);
        
        snapnow;
        close all;
        
        sprintf('Ca events, subtree %s is the only one that active', gRoi.Nodes.Name{classes(cl)})
        plotHMEvents(allEventsTable(eventsCl4Location, :), roiActivity_comb, selectedROI, [outputpath, '\only_',  gRoi.Nodes.Name{classes(cl)}], roiSortedByCluster, selectedROISplitDepth1, gRoi, roiTreeDistanceMatrix, index_apical);
        
        snapnow;
        close all;
        
        sprintf('Ca events, subtree %s is the only one that active', gRoi.Nodes.Name{classes(cl+1)})
        plotHMEvents(allEventsTable(eventsCl5Location, :), roiActivity_comb, selectedROI, [outputpath, '\only',  gRoi.Nodes.Name{classes(cl+1)}], roiSortedByCluster, selectedROISplitDepth1, gRoi, roiTreeDistanceMatrix, index_apical);
        
        snapnow;
        close all;
    end
    
    for cl = 1:length(classes)
        allEventsTable.(gRoi.Nodes.Name{classes(cl)}) = eventSplit(:, cl); 
    end
    
    writetable(allEventsTable, [outputpath, 'SummaryEventsWithP.csv']);
end

function plotHMEvents(allEventsTable, roiActivity_comb, selectedROI, outputpath, roiSortedByCluster, selectedROISplitDepth1, gRoi, roiTreeDistanceMatrix, index_apical)
    vectorToCompar = [];
    for ca_i = 1:size(allEventsTable, 1)
        vectorToCompar = [vectorToCompar, allEventsTable.start(ca_i):allEventsTable.pks(ca_i)];
    end
    
    activityMatrix = zeros(length(selectedROI), length(selectedROI), 2);
    for roi_i = 1:length(selectedROI)
        for roi_k = (roi_i):length(selectedROI)
            if (roi_k == roi_i)
                activityMatrix(roi_i, roi_i, 1) = 1;
            else
                currentROIActivity = roiActivity_comb(vectorToCompar, roi_i);
                secROIActivity = roiActivity_comb(vectorToCompar, roi_k);
                
                if isempty(vectorToCompar)
                    activityMatrix(roi_i, roi_k, 1) = nan;
                    activityMatrix(roi_k, roi_i, 1) = nan;
                else   
                    [corrEventsPeaks, pvalEventsPeaks] = corr([currentROIActivity, secROIActivity], 'type', 'Pearson');
                    activityMatrix(roi_i, roi_k, :) = [corrEventsPeaks(1, 2), pvalEventsPeaks(1,2)];
                    activityMatrix(roi_k, roi_i, :) = [corrEventsPeaks(1, 2), pvalEventsPeaks(1,2)];
                end         
            end
        end 
    end
    
    figDist = figure;
    hold on;
    title({'ROI Activity Distance', ['events:', num2str(size(allEventsTable, 1))]});
    xticks(1:length(selectedROI));
    yticks(1:length(selectedROI));
    aM = squeeze(activityMatrix(roiSortedByCluster, roiSortedByCluster, 1));
    m = imagesc(aM);
    colorbar
    cmap = jet();

    colormap(cmap);

    set(m,'AlphaData',~isnan(aM))

    for index_roi = 1:length(selectedROI)
        labelsNames(index_roi) = {sprintf('roi%d', sscanf(selectedROI{index_roi}, 'roi%d'))};
    end

    xticklabels(labelsNames(roiSortedByCluster));
    xtickangle(90);
    yticklabels(labelsNames(roiSortedByCluster));
    picNameFile = [outputpath, '\ActivityHMBYEventsSplit'];
    mysave(figDist, picNameFile);  

    plotROIDistMatrixTreeVSActivity(size(allEventsTable, 1), gRoi, outputpath, selectedROISplitDepth1, selectedROISplitDepth1, roiTreeDistanceMatrix, activityMatrix, true, 'Pearson', 'Split', selectedROI, index_apical, 'Correlation', '');
end

function plotHMEventsP(allEventsTable, roiActivity_comb, selectedROI, outputpath, roiSortedByCluster, selectedROISplitDepth1, gRoi, roiTreeDistanceMatrix, index_apical) 
    eventsMatrix = reshape(cell2mat(allEventsTable.roisEvent), 20, []);
    
    activityMatrix = zeros(length(selectedROI), length(selectedROI));
    for roi_i = 1:length(selectedROI)
        for roi_k = (roi_i):length(selectedROI)
%             if (roi_k == roi_i)
%                 activityMatrix(roi_i, roi_i) = 1;
%             else
                r_calc = sum(eventsMatrix(roi_i, :) & eventsMatrix(roi_k, :)) ;
                activityMatrix(roi_i, roi_k) = r_calc;
                activityMatrix(roi_k, roi_i) = r_calc;        
%             end
        end 
    end
    
    figDist = figure;
    hold on;
    title({'ROI Activity Distance', ['events:', num2str(size(allEventsTable, 1))]});
    xticks(1:length(selectedROI));
    yticks(1:length(selectedROI));
    m = imagesc(activityMatrix(roiSortedByCluster, roiSortedByCluster));
    colorbar
    cmap = jet();

    colormap(cmap);

    set(m,'AlphaData',~isnan(activityMatrix(roiSortedByCluster, roiSortedByCluster)))

    for index_roi = 1:length(selectedROI)
        labelsNames(index_roi) = {sprintf('roi%d', sscanf(selectedROI{index_roi}, 'roi%d'))};
    end

    xticklabels(labelsNames(roiSortedByCluster));
    xtickangle(90);
    yticklabels(labelsNames(roiSortedByCluster));
    picNameFile = [outputpath, '\ActivityHMBYEventsSplit'];
    mysave(figDist, picNameFile);  

    plotROIDistMatrixTreeVSActivity(size(allEventsTable, 1), gRoi, outputpath, selectedROISplitDepth1, selectedROISplitDepth1, roiTreeDistanceMatrix, activityMatrix, true, 'Pearson', 'Split', selectedROI, index_apical, 'Correlation', '');
end
  