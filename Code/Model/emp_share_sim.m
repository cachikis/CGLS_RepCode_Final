%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: emp_share_sim.m
% Author: Craig A. Chikis
% Date: 09/15/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = emp_share_sim(m) 

    panel_use = m.panel_save((m.panel_save.Time >= max(m.panel_save.Time) - 10 + 1), :); 

    lt10 = sum(panel_use.emp(round(panel_use.Age) <= 10));
    lt5 = sum(panel_use.emp(round(panel_use.Age) <= 5));
    lt1 = sum(panel_use.emp(round(panel_use.Age) <= 1));
    lt3 = sum(panel_use.emp(round(panel_use.Age) <= 3));  
    tot = sum(panel_use.emp); 

    m.share1 = lt1/tot;
    m.share3 = lt3/tot;
    m.share10 = lt10/tot;
    m.share5 = lt5/tot;


end 