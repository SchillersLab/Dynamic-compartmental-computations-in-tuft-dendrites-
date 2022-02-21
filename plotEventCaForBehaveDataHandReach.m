function plotEventCaForBehaveDataHandReach(behaveStartTime,tr_frame_count, allEventsTable, clusterCount, outputpath, behaveName, trials_label, split_trialsLabel)   
    for i = 0:clusterCount
        if i == 0
            titleText = {'All Ca Events', ['Aligned By: ', behaveName]};
            
            currentCaTable = allEventsTable(allEventsTable.clusterByH ~= 1, :);
        else
            currentCaTable = allEventsTable(allEventsTable.clusterByH == i,:);
            titleText = {['Cluster ', num2str(i),  ' Ca Events'], ['Aligned By: ', behaveName]};
        end
        
        behaveStartTimeFix = tr_frame_count * (currentCaTable.tr_index - 1) + behaveStartTime(currentCaTable.tr_index);
        eventsStart = currentCaTable.start - behaveStartTimeFix;
          
        eventsEnd = min(currentCaTable.event_end - behaveStartTimeFix, currentCaTable.pks - behaveStartTimeFix + 50);
        
        minEventStart = min(eventsStart);
        maxEventEnd = max(eventsEnd);
        
        timeEvents =  minEventStart:maxEventEnd; 
        
        ca_activitySummary = zeros(clusterCount, length(timeEvents));
        
        fig = figure;
        hold on;
             
        fig.Position = [fig.Position(1), fig.Position(2)-100 ,fig.Position(3) + 50,fig.Position(4) + 100];
        
        s1 = subplot(8, 1, 1:5);
        hold on;
        
        title(titleText); 
          
        leg = [];
        for k = 1:clusterCount
            leg(k) = plot(0,0, 'Color', getClusterColor(k), 'LineWidth', 1.5);
            legColor(k) = {['cluster ' num2str(k)]};
        end
       
        for ca_i = 1:size(currentCaTable, 1)
            if split_trialsLabel ~= 0 & trials_label(currentCaTable.tr_index(ca_i)) ~= split_trialsLabel
               continue;
            end
               
            if isempty(behaveStartTime(currentCaTable.tr_index(ca_i))) | behaveStartTime(currentCaTable.tr_index(ca_i)) == 0 
                continue;
            end
            
            plot(ones(1,9)*eventsStart(ca_i), ca_i:0.1:(ca_i + 0.8), 'Color', getClusterColor(currentCaTable.clusterByH(ca_i)), 'LineWidth', 2.5);
            
            timeLocation = find(timeEvents == eventsStart(ca_i),1);
            ca_activitySummary(currentCaTable.clusterByH(ca_i), timeLocation) = ca_activitySummary(currentCaTable.clusterByH(ca_i), timeLocation) + 1;
        end
        
        
%         xlim([max(-100,minEventStart), min(100,maxEventEnd)]);
        xlim([-100, 100]);
        
        ylim([0, size(currentCaTable, 1) + 1]);
        plot(zeros(1, size(currentCaTable, 1)), 1:size(currentCaTable, 1), '--k');
        plot(ones(1, size(currentCaTable, 1))*-3, 1:size(currentCaTable, 1), '--k');
        plot(ones(1, size(currentCaTable, 1))*3, 1:size(currentCaTable, 1), '--k');
        ylabel('Ca Event#');
        xlabel('Frame#');
        legend(leg, legColor);
        legend('Location', 'best')
        
%         legend('show');
        
        s2 = subplot(8, 1, 7:8);
        hold on;
        
        totalEventsSum = length(currentCaTable.clusterByH);
        
        totalEventsClusterSum(1) = sum(currentCaTable.clusterByH == 1);
        totalEventsClusterSum(2) = sum(currentCaTable.clusterByH == 2);
        totalEventsClusterSum(3) = sum(currentCaTable.clusterByH == 3);
        totalEventsClusterSum(4) = sum(currentCaTable.clusterByH == 4);
        
        if (i == 0)
            plot(timeEvents, sum(ca_activitySummary)./totalEventsSum, 'Color', 'k');
        
            for index_cl = 1:clusterCount
                plot(timeEvents, ca_activitySummary(index_cl, :)./totalEventsClusterSum(index_cl), 'Color', getClusterColor(index_cl));
            end
        else
            plot(timeEvents, (ca_activitySummary(i, :))./totalEventsClusterSum(i), 'Color', getClusterColor(i));
        end
        
        ylim([0,1]);
        
        limitsY = ylim;
        plot(zeros(1, 2), [0,1], '--k');
        plot([-3,-3], [0,1], '--k');
        plot([3,3], [0,1], '--k');
        
        ylabel('Ca Event# Summary');
        xlabel('Frame#');
%         xlim([max(-100,minEventStart), min(100,maxEventEnd)]);
        xlim([-100, 100]);
        
        linkaxes([s1, s2], 'x');
        
        mysave(fig, [outputpath, '\BehaveAlignedHandreach\cluster_', num2str(i), '_', behaveName]);
    end
end

function c_cluster = getClusterColor(clusterNum)
    switch clusterNum
        case 1
            c_cluster = [145, 30, 180] ./ 255;
        case 2
            c_cluster = [0, 130, 200] ./ 255;
        case 3
            c_cluster = [230, 25, 75] ./ 255;
        case 4
            c_cluster = [245, 130, 48] ./ 255;
    end
end