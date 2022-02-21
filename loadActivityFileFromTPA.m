function [roiActivity, roiActivityNames, tr_frame_count] = loadActivityFileFromTPA(activityfileTPAFolder, selectedROI, outputpath)
    TPAList = dir(strcat(activityfileTPAFolder,'\TPA*.mat'));

    fileNumRoi = length(TPAList);
     
    dataSize = 0;
    activity = table;
    for trialInd = 1:fileNumRoi
        usrData                    = load(fullfile(TPAList(trialInd).folder, TPAList(trialInd).name));
        if ~isfield(usrData, 'strROI')
            error(' has no ROI');
        end
        for m = 1:length(usrData.strROI)
            currentData = [];
            
            currentROIName = sprintf('roi%05d', extractROIstr(usrData.strROI{m}.Name));
            
            if usrData.strROI{m}.Name(end-1) == '_'
                currentROIName = sprintf('%s%s', currentROIName,  usrData.strROI{m}.Name(end));
            end
                       
            indexROI = find(strcmp(activity.Properties.VariableNames, currentROIName), 1);
            if isempty(indexROI)&& trialInd ~= 1
                    error('ROI Not exists in all trials');
            end
            
            % match new format to old format, the deltaF over F is saved in Data(:,2)
            % instead of procROI
            if ~isfield(usrData.strROI{m}, 'Data') && ~isprop(usrData.strROI{m}, 'Data')
                if ~isfield(usrData.strROI{m}, 'procROI')
                    error(' unfamiliar TPA file, cannot extract data');
                else
                    currentData = usrData.strROI{m}.procROI;
                end
            else
                currentData = usrData.strROI{m}.Data(:,2);
            end
            
            if dataSize == 0
                dataSize = length(currentData);
            elseif dataSize ~= size(currentData, 1)
                warning('ROI not the same data size!!!!!!!!!!!!!!!!');
                break;
            end
            
            activity.(currentROIName)(((trialInd - 1) * (length(currentData)) + 1):((trialInd) * (length(currentData)))) = currentData;        
        end
    end

    tr_frame_count = dataSize;
    
    roinames = activity.Properties.VariableNames;
    roi_index = 1;
    for index = 1:length(selectedROI)
        sel_results = [];
        if contains(selectedROI(index), '&')
            str_roi = strsplit(selectedROI{index}, '&');
            sel_results(1) = find(contains(roinames, str_roi(1)));
            sel_results(2) = find(contains(roinames, str_roi(2)));
        else
            sel_results = find(contains(roinames, selectedROI{index}));
        end
         
        if ~isempty(sel_results)
            if length(sel_results) == 2
                currentROIActivity = (activity.(roinames{sel_results(1)}) + activity.(roinames{sel_results(2)})) ./ 2; 
                currentName = selectedROI{index};
            elseif length(sel_results) == 1
                currentROIActivity = activity.(roinames{sel_results(1)});
                currentName = roinames{sel_results(1)};
            else
                error();
            end
            
            roiActivity(:, roi_index) = currentROIActivity;
            roiActivityNames{roi_index} = currentName;
            roi_index = roi_index + 1;
        end
        
%         if isempty(selectedROI) || ~isempty(find(strcmpi(selectedROI, roinames{index}), 1))
%            currentROIActivity = activity.(roinames{index});
%            
%            roiActivity(:, roi_index) = currentROIActivity;
%            roiActivityNames{roi_index} = roinames{index};
%            roi_index = roi_index + 1;
%         end
    end
    
%     writetable(activity,fullfile(outputpath, 'activityFileAsTable.csv'))
end