function GLM_SVM(selectedROI, cont, inds, roiNamesList, timeSeg, selectedROISplitDepth1, outputpath)
   corrContribution = [];
   selectedROISplitDepth1_sub = [];
    for i_roi = 1:length(selectedROI)        
        roiNumI=sscanf(selectedROI{i_roi}, 'roi%05d');
        locationI_glm = find(roiNumI == roiNamesList);            
        contLocationI = find(inds{timeSeg} == locationI_glm);
            
        if ~isempty(contLocationI)
            corrContribution(end+1, :) = cont{timeSeg}(contLocationI, :);
            selectedROISplitDepth1_sub(end+1) = selectedROISplitDepth1(i_roi);
        end                       
    end
    
    selectedROISplitDepth1_sub = selectedROISplitDepth1_sub';
    
    [T, ACC2D_depth1] = evalc("svmClassifyAndRand(corrContribution, selectedROISplitDepth1_sub, selectedROISplitDepth1_sub, 10, '', 1, 0)");  
    
    chanceCalc = hist(selectedROISplitDepth1_sub, unique(selectedROISplitDepth1_sub));
    chanceCalc = chanceCalc/sum(chanceCalc);
    ch_level = max(chanceCalc);
    
    save([outputpath, '\svm_glmVector.mat'], 'ACC2D_depth1', 'ch_level');
end
               