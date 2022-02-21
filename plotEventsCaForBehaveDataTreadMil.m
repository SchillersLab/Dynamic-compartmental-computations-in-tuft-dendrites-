function plotEventsCaForBehaveDataTreadMil(speedActivity, accelActivity, allEventsTable, clusterCount, outputpath, win_plot, BehaveDataTreadmil)
              
    for i = 0:clusterCount
        if i == 0
            currentCaTable = allEventsTable;
        else
            currentCaTable = allEventsTable(allEventsTable.clusterByH == i,:);
        end
        
        
        plotBehaveTreadMil(speedActivity, currentCaTable, win_plot, 'Velocity', i, outputpath)
        plotBehaveTreadMil(accelActivity, currentCaTable, win_plot, 'Acceleration', i, outputpath)
        
    end
end

function plotBehaveTreadMil(behave, currentCaTable, win_plot, behaveType, cluster_i, outputpath)
    tBehave = (1:length(behave)) - currentCaTable.start;
    
    minStart = min(tBehave(:, 1)); 
    maxEnd = max(tBehave(:, end));
    
    timeEvents =  minStart:maxEnd; 
    behave_Summary = zeros(size(currentCaTable, 1), length(timeEvents));
        
    for ca_i = 1:size(currentCaTable, 1)
        timeSLocation = find(timeEvents == tBehave(ca_i, 1),1);        
        timeELocation = find(timeEvents == tBehave(ca_i, end),1);
        
        behave_Summary(ca_i, timeSLocation:timeELocation) = behave;        
    end
    
    fig = figure;
    hold on;

    s1 = subplot(2, 1, 1);
    imagesc([timeEvents(1),timeEvents(end)], [1, size(currentCaTable, 1)], behave_Summary);
    cmap = jet();
    colormap(cmap);
    title({['Behave ' behaveType ' aligned To Ca Events'], ['cluster : ' num2str(cluster_i)]});
    xlabel('Frame');
    ylabel('Event#');
    hold on;
    plot(zeros(1, size(currentCaTable, 1)), 1:size(currentCaTable), '--k', 'LineWidth', 1.5);
    
    behave_Summary(isnan(behave_Summary)) = 0;
    
    if size(behave_Summary,1) == 1
        meanBehave = behave_Summary;
        stdBehave = behave_Summary;
    else
        meanBehave = mean(behave_Summary);
        stdBehave = std(behave_Summary);
    end
    
    totalSizeBefore = length(meanBehave);
    remBinSize = rem(length(meanBehave), 10);
    meanBehave = mean(reshape(meanBehave(1:(totalSizeBefore-remBinSize)), 10, []));
    stdBehave = mean(reshape(stdBehave(1:(totalSizeBefore-remBinSize)), 10, []));
    timeEvents = mean(reshape(timeEvents(1:(totalSizeBefore-remBinSize)), 10, []));
    
    s2 = subplot(2, 1, 2);
    hold on;
    H = mseb(timeEvents, meanBehave, stdBehave, [], 0.5);
    plot(zeros(1, 4), -0.05:0.05, '--k', 'LineWidth', 1.5);
    
    title({['Behave ' behaveType ' aligned To Ca Events'], ['cluster : ' num2str(cluster_i)]});
    xlabel('Frame');
    ylabel([behaveType, ' average']);
   
    mysave(fig, [outputpath, '\TreadMilAndCa\Behave_' behaveType '_AlignedToCaEvents_cluster' num2str(cluster_i)]);  
end
