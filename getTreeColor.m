function color = getTreeColor(methodT, treeNum, isMainSubTree, subTreeCountMain)
    if ~isMainSubTree
        treeNum = treeNum + subTreeCountMain;
    end
    
    switch methodT
        case 'ND'
            color = [0,0,0];
        case 'main'
            color = [0.82,0.22,0.31];
        case 'between'
            color = [25 110 180] ./ 255;
        case 'within'
            switch treeNum
                case 1
                    color = [1.00,0.60,0.20];
                case 2
                    color = [0.44,0.75,0.43];
                case 3
                    color = [0.68,0.44,0.71];
                case 4
                    color = [0.69,0.40,0.24];
                case 5
                    color = [0.87,0.81,0.43]; 
                case 6
                    color = [0.97,0.56,0.77];
                case 7
                    color = [176,196,222] ./ 255;
                case 8
                    color = [70, 130, 180] ./ 255;
                case 9
                    color = [115, 194, 251] ./ 255;
                case 10
                    color = [250, 128, 114] ./ 255;
                otherwise
                    color = [0,0,0];
            end
    end
end
