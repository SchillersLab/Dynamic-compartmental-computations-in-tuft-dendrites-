function [pVal, out] = calcMantelPerEvent(clusterType, curEventPearson, roiTreeDistanceMatrix, outputpathCurr, eventTable, clusterCount)
    for e_i = 1:size(curEventPearson, 1)
        curActivityP = squeeze(curEventPearson(e_i, :, :));
        curActivityP(isnan(curActivityP)) = 0;
        curActivityP(curActivityP < 0) = 0;
        
        if sum(curActivityP == 0, 'all') == (size(curActivityP, 1)*size(curActivityP, 1) - size(curActivityP, 1)) 
            out(e_i) = 0;
            pVal(e_i) = 1;
            continue;
        end
        
        if size(curActivityP, 1) <= 2
            out(e_i) = 0;
            pVal(e_i) = 1;
            continue;
        end
        
        [out(e_i), pVal(e_i)] = bramila_mantel(1- abs(curActivityP), roiTreeDistanceMatrix, 5000, 'pearson');
    end
    
    tr = table(out', pVal', eventTable.event_name, eventTable.clusterByH, eventTable.clusterByRoiPrecantage,...
        'VariableNames', {'Static', 'pVal', 'eventName', 'clusterByH', 'clusterByP'});
    mkdir(fullfile(outputpathCurr, clusterType));
    
    writetable(tr, fullfile(outputpathCurr, clusterType, 'mentalResults.csv'));
end