function saveNewTPAFile(selectedROI, roiActivity_comb, tr_frame_count, output, activityfileTPAFolder)
    TPAList = dir(strcat(activityfileTPAFolder,'\TPA*'));

    fileNumRoi = length(TPAList);
    
    for trialInd = 1:fileNumRoi
        usrData = load(fullfile(TPAList(trialInd).folder, TPAList(trialInd).name));
        
        for index_roi = 1:length(usrData.strROI)
            findROI = find(selectedROI, sprintf('roi%05d', extractROIstr(usrData.strROI(index_roi).Name)));
            usrData.strROI(index_roi).Data(:,2) = roiActivity_comb(((trialInd -1) * tr_frame_count + 1):((trialInd) * tr_frame_count), findROI);
        end
        
        
    end

end