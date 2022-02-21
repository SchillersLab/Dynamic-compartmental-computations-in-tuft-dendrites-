function TestBehaveEventsCorr(tableResults1)
    k = 1;
    for i = 1:size(tableResults1)
        load([tableResults1.RunLocation{i}, '\runParametes.mat'], 'globalParameters');
        if globalParameters.isHandreach
            neuronActivityPathTPA = fullfile(globalParameters.TPAFolder, globalParameters.AnimalName, globalParameters.DateAnimal);
            [~, ~, trials_label] = loadBDAFile(neuronActivityPathTPA, globalParameters.BehavioralSamplingRate, globalParameters.ImageSamplingRate, globalParameters.ImageSamplingRate*globalParameters.time_sec_of_trial, globalParameters.behavioralDelay, globalParameters.toneTime);
            eventsTable = readtable([tableResults1.RunLocation{i}, '\eventsCaSummary.csv']);
            trIndexAll = 1:length(trials_label);
            failCluster2and3(k) = sum(eventsTable.tr_index(eventsTable.clusterByH ~= 1& eventsTable.clusterByH ~= 4) == trIndexAll(trials_label == 2), 'all');
            sucCluster2and3(k) = sum(eventsTable.tr_index(eventsTable.clusterByH ~= 1 & eventsTable.clusterByH ~= 4) == trIndexAll(trials_label == 1), 'all');
            k = k+1;
            indexOfTable(k-1) = i;
        end
    end
end