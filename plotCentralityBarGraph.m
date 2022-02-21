function  plotCentralityBarGraph(centH, selectedROI, roiSortedByCluster, clusterCount, colorMatrix1, outputpath, roiLinkage)
    centrality_measures = fieldnames(centH);
        
    for i = 1:length(centrality_measures)
       resultsSummary = [];
       
       labels = cell(1, clusterCount + 1);
       for k = -1:clusterCount
           labels(k+2) = {sprintf('cluster%d', k)};
           resultsSummary(:, k+2) = centH(k+2).(centrality_measures{i});
       end
       
       labels(2) = {replace(labels{2}, 'cluster0', 'All')};
       
       fig = figure;
       hold on;
       sgtitle(sprintf('Centrality Summary - %s', centrality_measures{i}));
       
       fig.Color = [1,1,1];
       sb1 = subplot(6, 1, 5:6);

       d = dendrogram(roiLinkage, length(selectedROI), 'Labels', selectedROI, 'reorder', roiSortedByCluster, 'Orientation', 'bottom');    
       xlim([0.5,length(selectedROI)+0.5]);

       set(gca, 'color', 'none');
       axis off;
       set(d, 'Color', [0,0,0]);
       set(d, 'LineWidth', 1.5);

       sb2 = subplot(6, 1, 1:3);
       hold on;
     
       b = bar(resultsSummary(roiSortedByCluster, 3:end));
       
       for b_index = 1:length(b)
           b(b_index).FaceColor = 'flat';
           b(b_index).CData = colorMatrix1(roiSortedByCluster, :);
           b(b_index).FaceAlpha = 0.15 * b_index;
       end
       
%        xlabel('ROI"s');
       ax = gca;
       ax.Box = 'off';

       ylabel('Centrality');
       xtickangle(90);
       xticks(1:length(selectedROI));
       xlim([0.5, length(selectedROI)+0.5]);
       xticklabels(selectedROI(roiSortedByCluster));
       
       leg = [];
       colorU = unique(colorMatrix1, 'rows');
       index_leg = 1;
       for c_u = 1:size(colorU, 1)
           for k = 2:clusterCount+1
                leg(index_leg) = plot(0,0, 'Color', [colorU(c_u, :), 0.15 * k],  'LineWidth', 1.5);
                legColor(index_leg) = labels(k+1);
                index_leg = index_leg + 1;
           end
       end
       
       legend(leg, legColor);
       legend('Location', 'bestoutside');
       
       
       set(sb2, 'Position', [sb2.Position(1), sb2.Position(2) - 0.15, sb2.Position(3), sb2.Position(4)]);
       set(sb1, 'Position', [sb2.Position(1), sb2.Position(2) - sb1.Position(4), sb2.Position(3), sb1.Position(4)]);    
       
       mysave(fig, fullfile(outputpath, centrality_measures{i}, 'SummaryResults'))
    end   
end