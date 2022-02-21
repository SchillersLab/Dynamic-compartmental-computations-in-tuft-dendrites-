function plotROIDistMatrixTreeVSActivityBlack(event_count, gRoi, outputpath, firstBranchROI,mainTreeBranchROI, roiTreeDistanceMatrix, roiActivityDistanceMatrix, do_corrtest, roiActivityDistanceFunction, roiActivityPeakSize, selectedRoi, index_apical, distType, clusterType)          
   %     Plot ROI Activity VS Tree Distance
    fig = figure;
    hold on;
    subplot(9,1,2:8);
    sgtitle({'ROI Activity VS Tree Distance'});
    
    % Create ylabel
%     xlabel({'Calcium Event ' distType, ['(', num2str(event_count),')']});
    ylabel({'Calcium Event ' distType, ['(', num2str(event_count),')']});

    % Create xlabel
%     ylabel({'Dendritic distance'});
    xlabel({'Dendritic distance'});

    resultsForGLM = {};
    
    for index = 1: size(roiTreeDistanceMatrix, 2)
        for secIndex = (index + 1): size(roiTreeDistanceMatrix, 2)            
            hold on;
            scat = scatter(roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex, 1),20, 'filled', 'MarkerFaceColor', [0,0,0]);
            
            scat.set('UserData', [selectedRoi{index} 'x' selectedRoi{secIndex}])

            resultsForGLM(end + 1, :) = {([selectedRoi{index} 'x' selectedRoi{secIndex}]),...
                roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex, 1), roiActivityDistanceMatrix(index, secIndex, 2)};  
        end
    end
    
    mdAll = fitglm(cell2mat(resultsForGLM(:,2)), cell2mat(resultsForGLM(:,3)));
    yfitAll = predict(mdAll, cell2mat(resultsForGLM(:,2)));
    plot(cell2mat(resultsForGLM(:,2)), yfitAll, 'Color', [119, 123, 126] ./ 255, 'LineStyle', '--');

    txtSubT = sprintf(' {\\color[rgb]{%f,%f,%f}R^2=%0.4f}\n\r',0,...
        0, 0, mdAll.Rsquared.Ordinary);

    pV = foldsCalculator(cell2mat(resultsForGLM(:,2)), cell2mat(resultsForGLM(:,2)), mdAll.Rsquared.Ordinary);
    
    txtSubT = strcat(txtSubT, sprintf(' pVal = %f', pV));
    title(txtSubT);
   
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn})
    
    fig.Position = [fig.Position(1)-250, fig.Position(2)-100, fig.Position(3) + 250, fig.Position(4) + 100];
     
    fileName1 = [outputpath, '\ActivityDistVSDendriticDistForROI_' roiActivityDistanceFunction ,'_eventsSize', roiActivityPeakSize, '_numofTreeDepth_black'];
    mysave(fig, fileName1);
 
    snapnow;
    fclose('all');
    close all;
end