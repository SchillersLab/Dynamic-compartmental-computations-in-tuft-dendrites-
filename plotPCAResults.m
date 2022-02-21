function plotPCAResults(classes, resultsData, roiActivityNames, colorMatrix1, colorMatrix2, outputpath, clusterType, classesD1, classesD2, gRoi, selectedROISplitDepth1, selectedROISplitDepth3)
    
    for i = 1:length(classes)
        nameC = ['cluster_', num2str(classes(i))];
        tmp = resultsData.(nameC);
        tmp(isnan(tmp(:, 1)), :) = [];
        
        [pcares.embedding, ~, vals] = pca(tmp);
        pcares.effectiveDim = max(getEffectiveDim(vals, 0.95), 3);
        pcares.eigs = vals;
        embedding = pcares.embedding(:,1:2);
        
        [T, ACC2D_depth1] = evalc("svmClassifyAndRand(embedding, selectedROISplitDepth1, selectedROISplitDepth1, 10, '', 1, 0)");
        [T, ACC2D_depth2] = evalc("svmClassifyAndRand(embedding, selectedROISplitDepth3, selectedROISplitDepth3, 10, '', 1, 0)");
        
        chanceCalc = hist(selectedROISplitDepth1, unique(selectedROISplitDepth1));
        chanceCalc = chanceCalc/sum(chanceCalc);
        
        figPCA = figure;
        hold on;
        
        leg = [];
        for k = 1:length(classesD1)
            leg(k) = plot(0,0, 'Color', getTreeColor('within',k, true), 'LineWidth', 1.5);
            legColor(k) = gRoi.Nodes.Name(classesD1(k));
        end
       
        
        for c = 1:length(roiActivityNames)
            scatter(pcares.embedding(c, 1), pcares.embedding(c,2), 'filled', 'MarkerEdgeColor', colorMatrix1(c, :), 'MarkerFaceColor', colorMatrix1(c, :), 'SizeData', 20);
        end
        
        xlabel('PC1');
        ylabel('PC2');
        title({sprintf('PCA %s, Type %s',nameC, clusterType), 'Depth 1',  sprintf('accuracy mean %f, std %f, chance %f,', ACC2D_depth1.mean, ACC2D_depth1.std, max(chanceCalc))});
        legend(leg, legColor);
        legend('Location', 'best')
        
        mysave(figPCA, [outputpath, '\PcaResults\depth1_' nameC '_', clusterType]);     
        
        figPCA = figure;
        hold on;
        
        leg = [];
        for k = 1:length(classesD2)
            leg(k) = plot(0,0, 'Color', getTreeColor('within',k, false, length(classesD1)), 'LineWidth', 1.5);
            legColor(k) = gRoi.Nodes.Name(classesD2(k));
        end
        
        for c = 1:length(roiActivityNames)
            scatter(pcares.embedding(c, 1), pcares.embedding(c,2), 'filled', 'MarkerEdgeColor', colorMatrix2(c, :), 'MarkerFaceColor', colorMatrix2(c, :), 'SizeData', 20);
        end
        
        xlabel('PC1');
        ylabel('PC2');
        legend(leg, legColor);
        legend('Location', 'best')
        
        title({sprintf('PCA %s, Type %s',nameC, clusterType), 'Depth 2', sprintf('accuracy mean %f, std %f', ACC2D_depth2.mean, ACC2D_depth2.std)});
     
        mysave(figPCA, [outputpath, '\PcaResults\depth2_' nameC '_', clusterType]);
    end   
end