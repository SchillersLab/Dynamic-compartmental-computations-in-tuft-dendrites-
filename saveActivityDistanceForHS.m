function saveActivityDistanceForHS(roiActivityDistanceMatrixByH, outputpathCurr, roiActivityNames, colorMatrix1, colorMatrix2)      
     save([outputpathCurr, '\ActivityMatrixForHS.mat'], 'roiActivityDistanceMatrixByH', 'roiActivityNames', 'colorMatrix1', 'colorMatrix2');
end       