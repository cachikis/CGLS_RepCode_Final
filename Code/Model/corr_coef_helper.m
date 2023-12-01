%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: corr_coef_helper.m
% Author: Craig A. Chikis
% Date: 10/06/2023
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = corr_coef_helper(x, y)

    out = corrcoef(x, y); 

    if (prod(size(out)) == 4)
        out = out(1, 2); 
    else
        error("corr_coef_helper: corrcoef() returned a matrix that is not 2x2")
    end

end