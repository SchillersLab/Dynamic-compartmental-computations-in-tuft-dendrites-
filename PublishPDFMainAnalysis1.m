function PublishPDFMainAnalysis1(globalParameters)
    mainRunnerNeuronTreeAndActivityAnalysis_V3(globalParameters);
    
    if globalParameters.isHandreach && all(strcmp(globalParameters.runByEvent, 'non'))
        RunnerHadasCode(globalParameters);
    end
end