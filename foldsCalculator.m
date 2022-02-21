function [p_value,r_2]  = foldsCalculator(x,y, trueR)
    r_2 = [];
    for i = 1:1000
        p = randperm(length(x));
        mdl = fitglm(x(p), y);
        r_2(end+1) = mdl.Rsquared.Ordinary;
    end
    
    p_value = sum(r_2 >= trueR)/ 1000;
end