%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: innovation_output.m
% Author: Craig A. Chikis
% Date: 08/27/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = innovation_output(m)


    m.uncond_IO = [nanmean(m.panel_save.theta), nanstd(m.panel_save.theta), prctile(m.panel_save.theta, [1,5,10,25,50,75,90,95,99])]; 




end