%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: emp_growth_dispersion.m
% Author: Craig A. Chikis
% Date: 09/03/2023
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = emp_growth_dispersion(m)
    
    m.emp_growth_dispersion = nanstd(m.panel_save.empgrowth); 
    m.emp_growth_gross_reallocation = nanmean(abs(m.panel_save.empgrowth)); 



    
end