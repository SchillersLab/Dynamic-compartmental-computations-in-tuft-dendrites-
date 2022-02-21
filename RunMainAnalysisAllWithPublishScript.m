addpath('C:\Users\Jackie\Desktop\HadasCode\AnalysisProject\MatlabAnalysis');

globalParameters.TPAFolder = 'H:\Dropbox (Technion Dropbox)\Yara\Layer 5_Analysis\';
globalParameters.MainFolder = '\\jackie-analysis\e\';
globalParameters.AnimalName = 'Flp2';
globalParameters.DateAnimal = '09.27.20_LowTuft';
globalParameters.swcFile = 'swcFiles_LowTuft_Flp2\neuron_1.swc';
globalParameters.neuronNumberName = 'N1';
globalParameters.RunnerDate = 'final';
globalParameters.RunnerNumber = 'Run1';
globalParameters.behaveType = 'no_behave';
globalParameters.analysisType = 'Pearson\SP';
globalParameters.treadmilFile = '';
globalParameters.trajFolderName = '2020-09-27_FLP2_LOWTUFT';
globalParameters.videoPath = '';

globalParameters.hyperbolicDistMatrixLocation = "";
% globalParameters.hyperbolicDistMatrixLocation = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName,...
%     globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
%     globalParameters.RunnerDate, globalParameters.RunnerNumber, 'HS_create\StructureTreeHS\matlab_matrixbernoulli_100_3000.mat'); 

globalParameters.glmResults = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName,...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'glmAnalysis', 'glmResultsSummary.mat');

% -------------------------------------------------
globalParameters.runMantelTestPerDistanceThreshold = false;
globalParameters.runMantelTestPerDistanceThreshold_only = false;
globalParameters.MantelTJump = 3;
globalParameters.MantelTJumpFolder = fullfile(globalParameters.MainFolder, 'Shay',...
    'MantelThresholdTestSummary', globalParameters.AnimalName, globalParameters.DateAnimal,...
    globalParameters.neuronNumberName, globalParameters.RunnerDate,globalParameters.RunnerNumber, globalParameters.behaveType, globalParameters.analysisType);

% ----------------------------------------------------

%     Can be Euclidean OR ShortestPath OR ShortestPathCost OR Branch OR HyperbolicDist_L OR HyperbolicDist_P  = ( Between 2 roi according to the tree path )
globalParameters.roiTreeDistanceFunction = 'ShortestPath';
globalParameters.costSP = 0;

globalParameters.DistType = 'Pearson';

globalParameters.costCluster = -1;

globalParameters.std_treadMilThreshold = 0.0002;
globalParameters.aftertonetime = 16;

globalParameters.doExtraAnalysis = true;

globalParameters.treadMilExtraPlot = false;

globalParameters.reRunClusterData = false;
   
globalParameters.behaveFrameRateTM = 100;

globalParameters.ImageSamplingRate = 30;
globalParameters.time_sec_of_trial = 12;
globalParameters.trialNumber = [1, 2];
globalParameters.BehavioralSamplingRate = 200;
globalParameters.behavioralDelay = 0;
globalParameters.toneTime = 4;

%     No events put non
globalParameters.runByEvent = {'non'};
globalParameters.isHandreach = true;
globalParameters.EventTiming = 'start';
globalParameters.runByNegEvent = false;
globalParameters.FirstEventAfter = {'non'};

%     FOR Hand Reach 
%     NO labels 0, 1 suc , 2 fail
globalParameters.split_trialsLabel = 0;
globalParameters.runBehaveLag = [-inf, inf];
globalParameters.do_eventsBetween = false;
globalParameters.doBehaveAlignedPlot = false;

globalParameters.excludeTrailsByEventCount.Name = 'non';
globalParameters.excludeTrailsByEventCount.countRange = [-inf, inf];

%     Can be WindowEventFULLPearson OR
%     WindoEventToPeakPearson OR PeaksPearson OR glm
%     OR WindoEventToPeakCov OR WindoEventPearsonPerEvent OR HS_Activity
globalParameters.roiActivityDistanceFunction = 'WindoEventToPeakPearson';

globalParameters.clusterCount = 4;
globalParameters.eventWin = 10;

% 0.01 
globalParameters.mean_aV = 0.01; 
globalParameters.aVForAll = 0.3;
globalParameters.aVFix.location = [];
globalParameters.aVFix.values = [];

globalParameters.sigmaChangeValue = 0;
    
globalParameters.sigmaFix.location = [];
globalParameters.sigmaFix.values = [];

% most of the time 1.5 is good
globalParameters.thresholdGnValue = 4;
globalParameters.thresholdGnFix.location = [];
globalParameters.thresholdGnFix.values = [];

globalParameters.runMLS = true;
globalParameters.isSimData = false;

globalParameters.excludeRoi = [];

globalParameters.apical_roi = [];

globalParameters.firstDepthCompare = 1;
globalParameters.secDepthCompare = -1;

globalParameters.outputpath = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, ...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
    globalParameters.RunnerDate,globalParameters.RunnerNumber, globalParameters.behaveType, globalParameters.analysisType);

globalParameters.outputpath = char(globalParameters.outputpath);
mkdir(globalParameters.outputpath);

warning('off','all');

codePar = 'PublishPDFMainAnalysis1(globalParameters)';

% publish('PublishPDFMainAnalysis1.m', 'showCode', false,...
%     'format','pdf', 'catchError', false, 'outputDir', globalParameters.outputpath, 'codeToEvaluate', codePar);

% PublishPDFMainAnalysis1(globalParameters);

% mainRunnerNeuronTreeAndActivityAnalysis_V3(globalParameters);

RunnerHadasCode(globalParameters);


clear;