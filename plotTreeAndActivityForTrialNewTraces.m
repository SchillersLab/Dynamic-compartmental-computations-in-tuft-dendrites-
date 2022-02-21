function plotTreeAndActivityForTrialNewTraces()
   eventsListTableLocation = '\\jackie-analysis\e\Shay\PN11\02.03.22_ETL_Treadmill\Analysis\N1\Structural_VS_Functional\final\Run1\no_behave\Pearson\SP\';
   outputpath = '\\jackie-analysis\e\Shay\PN11\02.03.22_ETL_Treadmill\Analysis\N1\Structural_VS_Functional\final\Run1\tracesSoma\';
   neuronTreePathSWC = '\\jackie-analysis\e\Shay\PN11\02.03.22_ETL_Treadmill\swcFiles\neuron_1.swc';
   mkdir(outputpath);
   eventsList = readtable([eventsListTableLocation,'\eventsCaSummary.csv']);
   ImagingSamplineRate = 10.4167;
   
   trialTime = 12;
   eventWindow = [30,45];
   roiActivityLocation = eventsListTableLocation;
   load([roiActivityLocation, '\roiActivityRawData.mat'], 'roiActivity', 'roiActivityNames');
   
   excludeRoi = [10];
  
   trialNumber = [9,10,11,14,42,46,48];
   EventNumber = [1,1,1,1,1,1,1];
 
   colorMapLim = [0,0.03];
   
   windowIndexArray = [];
   eventLocationIndex = [];
   totalTrialTime = ImagingSamplineRate*trialTime;
   
   t = linspace(0, 12 * (round(size(roiActivity, 1) ./ totalTrialTime)), size(roiActivity, 1));
   
   
    [gRoi, rootNodeID, selectedROITable] = loadSwcFile(neuronTreePathSWC, outputpath, false);
    
    for i = 1:length(excludeRoi)
        ex_results = contains(selectedROITable.Name, sprintf('roi%05d', excludeRoi(i)));
        
        if sum(ex_results) == 1
            selectedROITable(ex_results, :) = [];
        end
    end
    
    selectedROI = selectedROITable.Name;
    selectedROISplitDepth1 = ones(length(selectedROI), 1) * -1;
    selectedROISplitDepth1 = getSelectedROISplitBranchID(gRoi, 1, selectedROISplitDepth1, selectedROI, rootNodeID);   
  
    [~, roiSortedByCluster, roiLinkage] = calcROIDistanceInTree_ShortestPath(gRoi, selectedROITable, outputpath, selectedROISplitDepth1, false);
      
   
   for trialIndex = 1:length(trialNumber)
        currentTrial = trialNumber(trialIndex);
        
        pksLocation = eventsList.pks(eventsList.tr_index == currentTrial);
        pksLocation = pksLocation(EventNumber(trialIndex));
        PlotTreeAndActivityForCurrentTrial(colorMapLim, currentTrial, totalTrialTime, roiSortedByCluster, roiActivity, roiActivityNames, selectedROI, outputpath, pksLocation, eventWindow, t, roiLinkage);
   
        windowIndexArray = [windowIndexArray, ((currentTrial - 1) * totalTrialTime + 1) : (currentTrial * totalTrialTime)]; 
       
        for curIndex = 1:length(pksLocation)
            eventLocationIndex = [eventLocationIndex,...
                    max(pksLocation(curIndex)-eventWindow(1), 1): min(pksLocation(curIndex) + eventWindow(2), size(roiActivity,1))];
        end
   end
   
   windowIndexArray = floor(windowIndexArray);
     
   for index = 1:length(selectedROI)
        activitylocation = strcmpi(roiActivityNames, selectedROI{index});
        SelectedROIActivity(index, :) = roiActivity(windowIndexArray, activitylocation);
                
        SelectedROIActivityWindowEvent(index, :) = roiActivity(eventLocationIndex, activitylocation);
   end
    
    fig = figure;
    hold on;

    sb1 = subplot(1, 6,1);

    dendrogram(roiLinkage, 'Labels', selectedROI, 'reorder', roiSortedByCluster, 'Orientation', 'left');
    ylim([0.5, length(selectedROI) + 0.5]);
    set(gca, 'color', 'none');
    axis off;
    
    sb2 = subplot(1, 6, 2:5);
 
    revA = roiSortedByCluster(:);
    
    for j = 1:length(revA)
        plot(SelectedROIActivity(revA(j), :)+(j)*2, 'k');;hold on;
    end
    
%     imagesc(t(windowIndexArray), 1:length(revA) ,SelectedROIActivity(revA, :));
    ax = gca;
    ax.YAxisLocation = 'right'; 
    yticks(2:2:length(revA)*2);
    yticklabels(selectedROI(revA));
    colormap('jet');
   
    ylim([1.5, length(selectedROI)*2 + 4]);
    caxis(colorMapLim);
    set(sb2, 'Position', [sb1.Position(1) + sb1.Position(3), sb1.Position(2), sb2.Position(3), sb2.Position(4)]);
    set(sb1, 'Position', [sb1.Position(1) , sb1.Position(2)-0.01, sb1.Position(3), sb1.Position(4)-0.05]);
    % Create title
    title(['Trial ', num2str(trialNumber)]);
    
    mysave(fig, [outputpath, '\DendrogramROIShortestPathDistAndActivityROIFORTrial_', num2str(trialNumber), 'AllTrail']);
    
%     -------------------------------------------------------
    fig2 = figure;
    hold on;

    sb1 = subplot(1, 6,1);

    dendrogram(roiLinkage, 'Labels', selectedROI, 'reorder', roiSortedByCluster, 'Orientation', 'left');
    ylim([0.5, length(selectedROI) + 0.5]);
    set(gca, 'color', 'none');
    axis off;
    
    sb2 = subplot(1, 6, 2:5);
     
    for j = 1:length(revA)
        plot(SelectedROIActivityWindowEvent(revA(j), :)+(j)*2, 'k');;hold on;
    end
  
%     imagesc(t(eventLocationIndex), 1:length(revA), SelectedROIActivityWindowEvent(revA, :));
    ax = gca;
    ax.YAxisLocation = 'right'; 
    yticks(2:2:length(revA)*2);
    yticklabels(selectedROI(revA));
    colormap('jet');
    
    ylim([1.5, length(selectedROI)*2 + 4]);
    caxis(colorMapLim);
    set(sb2, 'Position', [sb1.Position(1) + sb1.Position(3), sb1.Position(2), sb2.Position(3), sb2.Position(4)]);
    set(sb1, 'Position', [sb1.Position(1) , sb1.Position(2)-0.01, sb1.Position(3), sb1.Position(4)-0.05]);
      
    % Create title
    title(['Trial ', num2str(trialNumber), ' Events']);
    
    mysave(fig2, [outputpath, '\DendrogramROIShortestPathDistAndActivityROIFORTrial_', num2str(trialNumber), 'OnlyEventInTrail']);  
end

function PlotTreeAndActivityForCurrentTrial(colorMapLim, trialNumber, totalTrialTime, roiSortedByCluster, roiActivity, roiActivityNames, selectedROI, outputpath, pksLocation, eventWindow, t, roiLinkage)
   windowIndexArray = floor(((trialNumber - 1) * totalTrialTime + 1) : (trialNumber * totalTrialTime)); 
    
   eventLocationIndex = [];
   for curIndex = 1:length(pksLocation)
            eventLocationIndex = [eventLocationIndex,...
            max(pksLocation(curIndex) - eventWindow(1), 1) : min(pksLocation(curIndex) + eventWindow(2), size(roiActivity, 1))]; 
   end
   
   for index = 1:length(selectedROI)
        activitylocation = strcmpi(roiActivityNames, selectedROI{index});
        SelectedROIActivity(index, :) = roiActivity(windowIndexArray, activitylocation);
                
        SelectedROIActivityWindowEvent(index, :) = roiActivity(eventLocationIndex, activitylocation);
    end
    
    fig = figure;
    hold on;

    sb1 = subplot(1, 6,1);

    dendrogram(roiLinkage, 'Labels', selectedROI, 'reorder', roiSortedByCluster, 'Orientation', 'left');
    ylim([0.5, length(selectedROI) + 0.5]);
    set(gca, 'color', 'none');
    axis off;
    
    sb2 = subplot(1, 6, 2:5);
 
    revA = roiSortedByCluster(:);
         
    for j = 1:length(revA)
        plot(SelectedROIActivity(revA(j), :)+(j)*2, 'k');hold on;
    end
%     imagesc(t(windowIndexArray), 1:length(revA), SelectedROIActivity(revA, :));
    ax = gca;
    ax.YAxisLocation = 'right'; 
    yticks(2:2:length(revA)*2);
    yticklabels(selectedROI(revA));
    colormap('jet');
    ylim([1.5, length(selectedROI)*2 + 4]);
    caxis(colorMapLim);
    set(sb2, 'Position', [sb1.Position(1) + sb1.Position(3), sb1.Position(2), sb2.Position(3), sb2.Position(4)]);
    set(sb1, 'Position', [sb1.Position(1) , sb1.Position(2)-0.01, sb1.Position(3), sb1.Position(4)-0.05]);
    
    % Create title
    title(['Trial ', num2str(trialNumber)]);
    
    mysave(fig, [outputpath, '\DendrogramROIShortestPathDistAndActivityROIFORTrial_', num2str(trialNumber), 'AllTrail']);
    
%     -------------------------------------------------------
    fig2 = figure;
    hold on;

    sb1 = subplot(1, 6,1);

    dendrogram(roiLinkage, 'Labels', selectedROI, 'reorder', roiSortedByCluster, 'Orientation', 'left');
    ylim([0.5, length(selectedROI) + 0.5]);
    set(gca, 'color', 'none');
    axis off;
    
    sb2 = subplot(1, 6, 2:5);
    
    for j = 1:length(revA)
        plot(SelectedROIActivityWindowEvent(revA(j), :)+(j)*2, 'k');hold on;
    end
%     imagesc(t(eventLocationIndex), 1:length(revA), SelectedROIActivityWindowEvent(revA, :));
    ax = gca;
    ax.YAxisLocation = 'right'; 
    yticks(2:2:length(revA)*2);
    yticklabels(selectedROI(revA));
    colormap('jet');
  
    ylim([1.5, length(selectedROI)*2 + 4]);
    caxis(colorMapLim);
    set(sb2, 'Position', [sb1.Position(1) + sb1.Position(3), sb1.Position(2), sb2.Position(3), sb2.Position(4)]);
    set(sb1, 'Position', [sb1.Position(1) , sb1.Position(2)-0.01, sb1.Position(3), sb1.Position(4)-0.05]);
      
    % Create title
    title(['Trial ', num2str(trialNumber), ' Events']);
    
    mysave(fig2, [outputpath, '\DendrogramROIShortestPathDistAndActivityROIFORTrial_', num2str(trialNumber), 'OnlyEventInTrail']);  
end
