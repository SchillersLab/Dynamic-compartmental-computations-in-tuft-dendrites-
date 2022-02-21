function [roiActivity, roiActivityNames] = loadActivityFile(activityfile, selectedROI)
    activity = readtable(activityfile);
    roinames = activity.Properties.VariableNames;
    roi_index = 1;
    for index = 1:length(roinames)
        if isempty(selectedROI) || ~isempty(find(strcmpi(selectedROI, roinames{index}), 1))
           currentROIActivity = activity.(roinames{index});
           
           roiActivity(:, roi_index) = currentROIActivity;
           roiActivityNames{roi_index} = roinames{index};
           roi_index = roi_index + 1;
           
%            figROI = figure;
%           
%            % Create axes
%             axes1 = axes('Parent',figROI);
%             hold(axes1,'on');
% 
%             % Create ylabel
%             ylabel({'Activity'});
% 
%             % Create xlabel
%             xlabel({'Time'});
% 
%             % Create title
%             title({roinames{index}, 'Activity'});
% 
%             box(axes1,'on');
%             grid(axes1,'on');
% 
%             plot(1:length(currentROIActivity), currentROIActivity);
%             
        end
    end
end