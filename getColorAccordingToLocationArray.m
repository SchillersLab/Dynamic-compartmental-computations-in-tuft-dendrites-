function color = getColorAccordingToLocationArray(location, secLocation)
    if (location == secLocation)
        switch(location)
            case 1
                color = [233, 150, 122] ./ 255;
            case 2
                color = [143, 188, 143] ./ 255;
            case 3
                color = [72, 61, 139] ./ 255;
            case 4
                color = [0, 191, 255] ./ 255;
            case 5 
                color = [255, 105, 180] ./ 255;
            case 6
                color = [32, 178, 170] ./ 255;
            case 7
                color = [244, 164, 96] ./ 255;
            case 8
                color = [199, 21, 133] ./ 255;
            case 9
                color = [70, 130, 180] ./ 255;
            case 10
                color = [0, 0, 255] ./ 255;
            case 11
                color = [255, 255, 0] ./ 255;
            case 12
                color = [0, 255, 0] ./ 255;
        end
    else
        color = [rand(1), 0, 0];
    end
end