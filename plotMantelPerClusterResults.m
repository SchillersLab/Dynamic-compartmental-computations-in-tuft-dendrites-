function plotMantelPerClusterResults(outputpath, outM, pValM, clusterCount, otherPath)
    labels(1) = {'All'};
    
    for i = 1:clusterCount
        labels(i+1) = {sprintf('cluster%d', i)};        
    end
    
    fig = figure;
    sgtitle('Mental Test Per Cluster');
    hold on;
    subplot(8,1, 2:8);
    
    Y = [outM(2:end)]';
    b = bar(Y);
    ylim([-0.1, 1]);
    xticks(1:(clusterCount + 1));
    xticklabels(labels);
    xtickangle(90);
    xlabel('cluste#');
    ylabel('Mantel Test ');
    
    colorB(1, :) = [137, 207, 239] ./ 255;
    colorB(2, :) = [77, 81, 109] ./ 255;
%     labelsB = labels;
    
    text(1:length(Y),Y,num2str(pValM(2:end)'),'vert','bottom','horiz','center'); 
    
    for i = 1:length(b)
        b(i).FaceColor = colorB(1, :);        
    end
    
    mysave(fig, fullfile(outputpath, 'PerClusterMantelPlot'));
    
    mysave(fig, fullfile(otherPath, 'PerClusterMantelPlot'));
    
    matV = [pValM(2:end)', outM(2:end)'];
    cellV = num2cell(matV);
    cellV(:, end+1) = labels;
    tr = cell2table(cellV);
    tr.Properties.VariableNames = {'PValue', 'MantelValue', 'Cluster'};
    tr.Properties.RowNames = labels;
    writetable(tr, fullfile(outputpath, 'PerClusterMantelTable.csv'));
    
    
    writetable(tr, fullfile(otherPath, 'PerClusterMantelTable.csv'));
end