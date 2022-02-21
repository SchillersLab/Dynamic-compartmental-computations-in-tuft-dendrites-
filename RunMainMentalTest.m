addpath('C:\Users\Jackie\Desktop\HadasCode\AnalysisProject\MatlabAnalysis');

globalParameters.MainFolder = '\\jackie-analysis\e\';
globalParameters.AnimalName = 'SM04';
globalParameters.DateAnimal = '11.24.19_Tuft';
globalParameters.swcFile = 'SM04_11.24.19_Tuft_swcFiles\neuron_2.swc';
globalParameters.neuronNumberName = 'N2';
globalParameters.RunnerDate = '2-11-20';
globalParameters.RunnerNumber = 'Run2';
globalParameters.behaveType = 'no_behave';
globalParameters.analysisType = 'Pearson\SP';

globalParameters.hyperbolicDistMatrixLocation = "";
% globalParameters.hyperbolicDistMatrixLocation = fullfile(globalParameters.MainFolder, 'Shay', globalParameters.AnimalName,...
%     globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
%     globalParameters.RunnerDate, globalParameters.RunnerNumber, 'HS_create\StructureTreeHS\matlab_matrixbernoulli_100_3000.mat'); 

%     Can be Euclidean OR ShortestPath OR HyperbolicDist_L OR HyperbolicDist_P  = ( Between 2 roi according to the tree path )
globalParameters.roiTreeDistanceFunction = 'ShortestPath';

globalParameters.isFlipMapStructure = true;
globalParameters.isFlipMapActivity = false;

globalParameters.clusterCount = 4;

globalParameters.excludeRoi = [35];

globalParameters.apical_roi = [];

globalParameters.firstDepthCompare = 1;
globalParameters.secDepthCompare = 2;
globalParameters.thDepthCompare = 3;

globalParameters.outputpath = fullfile(globalParameters.MainFolder, 'Shay' , globalParameters.AnimalName, ...
    globalParameters.DateAnimal, 'Analysis', globalParameters.neuronNumberName, 'Structural_VS_Functional',...
    globalParameters.RunnerDate,globalParameters.RunnerNumber, globalParameters.behaveType, globalParameters.analysisType);

globalParameters.outputpath = char(globalParameters.outputpath);
mkdir(globalParameters.outputpath);

warning('off','all');

mainRunnerMentalTest(globalParameters);

clear;