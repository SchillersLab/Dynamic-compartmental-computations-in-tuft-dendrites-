function createRoiEventsSpecificForOldAnalysis
    runnerPath = '\\jackie-analysis\e\Shay\4481\2019-03-07_4481_motor_aligned\Analysis\N3\Structural_VS_Functional\final\Run2\';
    behaveType = 'walkPosacceleration';

    load(fullfile(runnerPath, 'no_behave', 'Pearson', 'SP', 'roiActivityRawData.mat'), 'roiActivityNames');
    tableEventsBehave = readtable(fullfile(runnerPath, behaveType, 'Pearson', 'SP' ,'eventsCaSummary.csv'));
    
    location = contains(fieldnames(tableEventsBehave), 'roisEvent_');
    
    for e_i = 1:size(tableEventsBehave, 1)
        tableEventsBehave.roisEvent(e_i) = {table2array(tableEventsBehave(e_i, location))};
    end
    
    calcRoiEventsSpecific(tableEventsBehave, roiActivityNames, fullfile(runnerPath, behaveType, 'Pearson', 'SP'), {behaveType})
end