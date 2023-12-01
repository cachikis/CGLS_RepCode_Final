%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: sim_wrap.m
% Author: Craig A. Chikis
% Date: 12/20/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = sim_wrap(m)
    if (m.entrance)
        m = sim_entry(m); 
    else
        m = sim_noentry(m); 
    end
end