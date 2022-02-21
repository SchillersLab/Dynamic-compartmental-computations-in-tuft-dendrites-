function [pictureNames, resultsForGLM] = plotROIDistMatrixTreeVSActivity_Type2(event_count, gRoi, outputpath, firstBranchROI,mainTreeBranchROI, roiTreeDistanceMatrix, roiActivityDistanceMatrix, do_corrtest, roiActivityDistanceFunction, roiActivityPeakSize, selectedRoi, index_apical, distType, clusterType)
    classesM = unique(mainTreeBranchROI);
    classesF = unique(firstBranchROI);
    
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

    
%     leg = zeros(3, 1);
    leg = [];
%     legColor = []; 
    
    
    index_classes = 1;   
    classesColorName = {};
    resultsForGLM = {};
    
    for index = 1: size(roiTreeDistanceMatrix, 2)
        for secIndex = (index + 1): size(roiTreeDistanceMatrix, 2)
            
            if mainTreeBranchROI(index) ~= -1 && mainTreeBranchROI(index) == mainTreeBranchROI(secIndex)
                locationFirstGROI = strcmp(gRoi.Nodes.Name, selectedRoi(index));
                locationSecGROI = strcmp(gRoi.Nodes.Name, selectedRoi(secIndex));
                
                [p, ~] = shortestpath(gRoi,gRoi.Nodes.ID(locationFirstGROI),gRoi.Nodes.ID(locationSecGROI), 'Method', 'unweighted');
                
                if any(p == mainTreeBranchROI(index))
                    fFirstName = gRoi.Nodes(mainTreeBranchROI(index),:).Name{1};
                    color = getTreeColor('within', find(classesM == mainTreeBranchROI(index)), false, length(classesF));
                    colorName = {fFirstName};
                    
                    hold on;
            
            
                    if sum((index_apical == index) | (index_apical == secIndex)) > 0
        %                 scat = scatter(roiActivityDistanceMatrix(index, secIndex, 1), roiTreeDistanceMatrix(index, secIndex), '*', 'MarkerEdgeColor', color);
                        scat = scatter(roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex, 1), '*', 'MarkerEdgeColor', color);
                    else
        %                 scat = scatter(roiActivityDistanceMatrix(index, secIndex, 1), roiTreeDistanceMatrix(index, secIndex), 20, 'filled', 'MarkerFaceColor', color);
                        scat = scatter(roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex, 1), 20, 'filled', 'MarkerFaceColor', color);           
                    end

                    scat.set('UserData', [selectedRoi{index} 'x' selectedRoi{secIndex}])

                    colorName = colorName{1};
                    colorName = replace(colorName,{'&','-', '_'}, 'x'); 

                    if isempty(classesColorName) || ...
                            (sum(strcmp(classesColorName, colorName)) == 0)
                        classesColorName(index_classes) = {colorName};
                        classesColor(index_classes) = {color};

                        leg(index_classes) = plot(0,0, 'color', color, 'LineWidth', 2.5);
                        legColor(index_classes) = {colorName};

                        resultsT.(colorName) = [];
                        index_classes = index_classes + 1;
                    end

                    resultsForGLM(end + 1, :) = {([selectedRoi{index} 'x' selectedRoi{secIndex}]),...
                        roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex, 1), roiActivityDistanceMatrix(index, secIndex, 2), mainTreeBranchROI(index), ...
                        mainTreeBranchROI(secIndex), (mainTreeBranchROI(secIndex) == mainTreeBranchROI(index) && mainTreeBranchROI(index) ~= -1), roiActivityPeakSize};  

                    resultsT.(colorName)(end + 1, :) = [roiTreeDistanceMatrix(index, secIndex), roiActivityDistanceMatrix(index, secIndex, 1), roiActivityDistanceMatrix(index, secIndex, 2)];

                end            
               
            end
        end
    end
 
    if ~exist('resultsT', 'var')
        fclose('all');
        close all;
        return;
    end
    
    f_name = fieldnames(resultsT);
    indexT = 1;
    txtSubT = '';
    r2 = {};
    for t = 1:length(f_name)
        if isempty(resultsT.(f_name{t}))
            continue;
        end
        locationColor = strcmp(classesColorName, f_name{t});
        
        if size(resultsT.(f_name{t}), 1) > 1
            mdl = fitglm(resultsT.(f_name{t})(:,1), resultsT.(f_name{t})(:,2));

            yfit = predict(mdl,resultsT.(f_name{t})(:,1));
    %         plot(resultsT.(f_name{t})(:,2), yfit, 'Color', classesColor{locationColor}, 'LineStyle', '--');
            plot(resultsT.(f_name{t})(:,1), yfit, 'Color', classesColor{locationColor}, 'LineStyle', '--');

            r2(indexT, 1) = {mdl.Rsquared.Ordinary};
            r2(indexT, 2) = f_name(t);

            txtSubT = strcat(txtSubT, sprintf(' {\\color[rgb]{%f,%f,%f}R^2=%0.4f}\n\r',classesColor{locationColor}(1),...
                classesColor{locationColor}(2), classesColor{locationColor}(3), mdl.Rsquared.Ordinary));

            indexT = indexT + 1;
        end
        writetable(array2table(resultsT.(f_name{t})),fullfile(outputpath, ['AVSDForROI_' num2str(t) '_Depth', num2str(length(classesM)), '_', f_name{t}, '.csv']));
    end
    
    title(txtSubT);
    
    if ~isempty(r2)
        writetable(cell2table(r2),fullfile(outputpath, ['Rsquared_' num2str(t) '_Depth', num2str(length(classesM)), '_', f_name{t}, '.csv']));
    end
    
    dcm_obj = datacursormode(fig);
    set(dcm_obj,'UpdateFcn',{@myupdatefcn})
    
    legend(leg, legColor);
    legend('Location','bestoutside')
    fig.Position = [fig.Position(1)-250, fig.Position(2)-100, fig.Position(3) + 250, fig.Position(4) + 100];
    roiwithstar = '*';
    
    for in = 1:length(index_apical)
       roiwithstar = [roiwithstar ' ' selectedRoi{index_apical(in)}]; 
    end
    
    annotation('textbox', [0, 0.2, 0, 0], 'string', roiwithstar, 'FitBoxToText','on')
%     ylim([0,1]);
       
    fileName1 = [outputpath, '\ActivityDistVSDendriticDistForROI_' roiActivityDistanceFunction ,'_eventsSize', roiActivityPeakSize, '_numofTreeDepth', num2str(length(classesM))];
    mysave(fig, fileName1);
    
    plotAndCalcStatisticTest(classesColorName, classesColor, resultsT, leg, legColor, outputpath, classesM, clusterType, roiActivityPeakSize, distType);
 
% ------------------------------------------------------------------------------
    
    %     coutI = 1;
    pictureNames = {fileName1};
    
    snapnow;
    fclose('all');
    close all;
    clear resultsT;
end