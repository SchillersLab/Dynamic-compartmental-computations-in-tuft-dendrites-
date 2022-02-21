function saveGraphForHS(gRoi, rootNodeID, outputpath, color1, color2, selectedRoiT)
    matrixForHS = getMatrixForDistanceHS(gRoi, rootNodeID); 
    matrixForHS_2 = getMatrixForDistanceHS_NoFix(gRoi, rootNodeID); 
    matrixForHS_3 = getMatrixForDistanceHS_Type2(gRoi, rootNodeID); 
    
    matrixForSP = getMatrixForDistanceSP(gRoi); 
    depthToSave = ceil(log2(size(gRoi.Nodes, 1) + 1));
    depthToSaveReal = (log2(size(gRoi.Nodes, 1) + 1));
    baches = size(gRoi.Nodes, 1);
    points_name = char(gRoi.Nodes.Name);
    
    colorMatrix1 = zeros(size(gRoi.Nodes, 1), 3);
    colorMatrix1(selectedRoiT.ID, :) = color1;
    
    colorMatrix2 = zeros(size(gRoi.Nodes, 1), 3);
    colorMatrix2(selectedRoiT.ID, :) = color2;
    
    save([outputpath, '\GraphAsMatrix.mat'], 'matrixForHS', 'depthToSave', 'baches', 'depthToSaveReal', 'points_name', 'matrixForSP', 'colorMatrix1', 'colorMatrix2');
    save([outputpath, '\GraphAsMatrix2.mat'], 'matrixForHS_2', 'depthToSave', 'baches', 'depthToSaveReal', 'points_name', 'matrixForSP', 'colorMatrix1', 'colorMatrix2');   
    save([outputpath, '\GraphAsMatrix3.mat'], 'matrixForHS_3', 'depthToSave', 'baches', 'depthToSaveReal', 'points_name', 'matrixForSP', 'colorMatrix1', 'colorMatrix2');   

    [matrixForHS_4, tempSize, matrixForSP] = getMatrixForDistanceHS_Type3(gRoi, rootNodeID); 
    
    for i = 1:tempSize
        points_name(end+1, 1:4) = char('temp');
        colorMatrix1(end+1, :) = [0, 0, 0];
        colorMatrix2(end+1, :) = [0, 0, 0];
    end
    
    depthToSave = ceil(log2(size(gRoi.Nodes, 1) + 1 + tempSize));
    depthToSaveReal = (log2(size(gRoi.Nodes, 1) + 1 + tempSize));
    baches = size(gRoi.Nodes, 1) + tempSize;
    
    save([outputpath, '\GraphAsMatrix4.mat'], 'matrixForHS_4', 'depthToSave', 'baches', 'depthToSaveReal', 'points_name', 'matrixForSP', 'colorMatrix1', 'colorMatrix2');   

end