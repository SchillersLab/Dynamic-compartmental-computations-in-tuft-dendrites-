function plotTreeAndActivityForTrial(trialNumber, totalTrialTime, roiSortedByCluster, roiActivity, roiActivityNames, selectedROI, outputpath, pksLocation, eventWindow, roiLinkage)
   windowIndexArray = [];
   eventLocationIndex = [];
   
   t = linspace(0, 12 * (round(size(roiActivity, 1) ./ totalTrialTime)), size(roiActivity, 1));
   
   for trialIndex = 1:length(trialNumber)
        currentTrial = trialNumber(trialIndex);
        PlotTreeAndActivityForCurrentTrial(currentTrial, totalTrialTime, roiSortedByCluster, roiActivity, roiActivityNames, selectedROI, outputpath, pksLocation, eventWindow, t, roiLinkage);
   
        windowIndexArray = [windowIndexArray, ((currentTrial - 1) * totalTrialTime + 1) : (currentTrial * totalTrialTime)]; 
        eventLocation = find(pksLocation >= ((currentTrial - 1) * totalTrialTime + 1) & pksLocation <= (currentTrial * totalTrialTime));
        
        for curIndex = eventLocation
            eventLocationIndex = [eventLocationIndex,...
                    eventWindow(curIndex, 1) : eventWindow(curIndex, 2)];
        end
   end
   
     
   for index = 1:length(selectedROI)
        activitylocation = strcmpi(roiActivityNames, selectedROI{index});
        SelectedROIActivity(index, :) = roiActivity(windowIndexArray, activitylocation);
                
        SelectedROIActivityWindowEvent(index, :) = roiActivity(eventLocationIndex, activitylocation);
   end
    
    fig = figure;
    hold on;

    sb1 = subplot(1, 6, 1:2);

    dendrogram(roiLinkage, 'Labels', selectedROI, 'reorder', roiSortedByCluster, 'Orientation', 'left');
    ylim([0.5, length(selectedROI) + 0.5]);
    set(gca, 'color', 'none');
    axis off;
    
    sb2 = subplot(1, 6, 3:5);
 
    revA = roiSortedByCluster(length(roiSortedByCluster):-1:1);
    imagesc(t(windowIndexArray), 1:length(revA) ,SelectedROIActivity(revA, :));
    ax = gca;
    ax.YAxisLocation = 'right'; 
    yticks(1:length(revA));
    yticklabels(selectedROI(revA));
    colormap('jet');
    
    set(sb2, 'Position', [sb1.Position(1) + sb1.Position(3), sb1.Position(2), sb2.Position(3), sb2.Position(4)]);
    set(sb1, 'Position', [sb1.Position(1) , sb1.Position(2), sb1.Position(3), sb1.Position(4)]);
    
    % Create title
    title(['Trial ', num2str(trialNumber)]);
    
    mysave(fig, [outputpath, '\DendrogramROIShortestPathDistAndActivityROIFORTrial_', num2str(trialNumber), 'AllTrail']);
    
%     -------------------------------------------------------
    fig2 = figure;
    hold on;

    sb1 = subplot(1, 6, 1:2);

    dendrogram(roiLinkage, 'Labels', selectedROI, 'reorder', roiSortedByCluster, 'Orientation', 'left');
    ylim([0.5, length(selectedROI) + 0.5]);
    set(gca, 'color', 'none');
    axis off;
    
    sb2 = subplot(1, 6, 3:5);
 
    imagesc(t(eventLocationIndex), 1:length(revA), SelectedROIActivityWindowEvent(revA, :));
    ax = gca;
    ax.YAxisLocation = 'right'; 
    yticks(1:length(revA));
    yticklabels(selectedROI(revA));
    colormap('jet');
    
    set(sb2, 'Position', [sb1.Position(1) + sb1.Position(3), sb1.Position(2), sb2.Position(3), sb2.Position(4)]);
    set(sb1, 'Position', [sb1.Position(1) , sb1.Position(2), sb1.Position(3), sb1.Position(4)]);
    
    % Create title
    title(['Trial ', num2str(trialNumber), ' Events']);
    
    mysave(fig2, [outputpath, '\DendrogramROIShortestPathDistAndActivityROIFORTrial_', num2str(trialNumber), 'OnlyEventInTrail']);  
end

function PlotTreeAndActivityForCurrentTrial(trialNumber, totalTrialTime, roiSortedByCluster, roiActivity, roiActivityNames, selectedROI, outputpath, pksLocation, eventWindow, t, roiLinkage)
   windowIndexArray = ((trialNumber - 1) * totalTrialTime + 1) : (trialNumber * totalTrialTime); 
   eventLocation = find(pksLocation >= windowIndexArray(1) & pksLocation <= windowIndexArray(end));
    
   eventLocationIndex = [];
   for curIndex = eventLocation
    eventLocationIndex = [eventLocationIndex,...
            eventWindow(curIndex, 1) : eventWindow(curIndex, 2)]; 
   end
   
   for index = 1:length(selectedROI)
        activitylocation = strcmpi(roiActivityNames, selectedROI{index});
        SelectedROIActivity(index, :) = roiActivity(windowIndexArray, activitylocation);
                
        SelectedROIActivityWindowEvent(index, :) = roiActivity(eventLocationIndex, activitylocation);
    end
    
    fig = figure;
    hold on;

    sb1 = subplot(1, 6, 1:2);

    dendrogram(roiLinkage, 'Labels', selectedROI, 'reorder', roiSortedByCluster, 'Orientation', 'left');
    ylim([0.5, length(selectedROI) + 0.5]);
    set(gca, 'color', 'none');
    axis off;
    
    sb2 = subplot(1, 6, 3:5);
 
    revA = roiSortedByCluster(length(roiSortedByCluster):-1:1);
    imagesc(t(windowIndexArray), 1:length(revA), SelectedROIActivity(revA, :));
    ax = gca;
    ax.YAxisLocation = 'right'; 
    yticks(1:length(revA));
    yticklabels(selectedROI(revA));
    colormap('jet');
    
    set(sb2, 'Position', [sb1.Position(1) + sb1.Position(3), sb1.Position(2), sb2.Position(3), sb2.Position(4)]);
    set(sb1, 'Position', [sb1.Position(1) , sb1.Position(2), sb1.Position(3), sb1.Position(4)]);
    
    % Create title
    title(['Trial ', num2str(trialNumber)]);
    
    mysave(fig, [outputpath, '\DendrogramROIShortestPathDistAndActivityROIFORTrial_', num2str(trialNumber), 'AllTrail']);
    
%     -------------------------------------------------------
    fig2 = figure;
    hold on;

    sb1 = subplot(1, 6, 1:2);

    dendrogram(roiLinkage, 'Labels', selectedROI, 'reorder', roiSortedByCluster, 'Orientation', 'left');
    ylim([0.5, length(selectedROI) + 0.5]);
    set(gca, 'color', 'none');
    axis off;
    
    sb2 = subplot(1, 6, 3:5);
 
    imagesc(t(eventLocationIndex), 1:length(revA), SelectedROIActivityWindowEvent(revA, :));
    ax = gca;
    ax.YAxisLocation = 'right'; 
    yticks(1:length(revA));
    yticklabels(selectedROI(revA));
    colormap('jet');
    
    set(sb2, 'Position', [sb1.Position(1) + sb1.Position(3), sb1.Position(2), sb2.Position(3), sb2.Position(4)]);
    set(sb1, 'Position', [sb1.Position(1) , sb1.Position(2), sb1.Position(3), sb1.Position(4)]);
    
    % Create title
    title(['Trial ', num2str(trialNumber), ' Events']);
    
    mysave(fig2, [outputpath, '\DendrogramROIShortestPathDistAndActivityROIFORTrial_', num2str(trialNumber), 'OnlyEventInTrail']);  
end
