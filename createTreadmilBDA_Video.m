function createTreadmilBDA_Video(globalParameters)
    neuronActivityPathTPA = fullfile(globalParameters.TPAFolder, globalParameters.AnimalName, globalParameters.DateAnimal);
    
    behaveFileTreadMillPath = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, globalParameters.DateAnimal, globalParameters.treadmilFile);
    
    eventsDetectionFolder = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, ...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
    globalParameters.RunnerDate,globalParameters.RunnerNumber, 'EventsDetection');
    
    load([eventsDetectionFolder, '\BehaveTreadMilResults'], 'BehaveDataTreadmil');
       
    saveDataFotGLMTreadMil_Type2(BehaveDataTreadmil, roiActivity, globalParameters.ImageSamplingRate, neuronActivityPathTPA, globalParameters.time_sec_of_trial,behaveFileTreadMillPath, globalParameters.videoPath);            

end