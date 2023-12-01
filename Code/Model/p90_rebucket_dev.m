%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: p90_rebucket_inp.m
% Author: Craig A. Chikis
% Date: 10/18/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [citation_track2, quality] = p90_rebucket_dev(citation_track2, uniquekeys, prctile_inp, method_prc)

    if (nargin < 4)
        method_prc = "exact";
    end

    func = @(x) prctile(x, prctile_inp, 'method', method_prc);
   


    citation_track2 = innerjoin(citation_track2, uniquekeys, 'keys', ["Time", "CPC", "Firm"]); 

    [~, uniquerows] = unique(citation_track2(:, ["wku", "CPC", "Time"]), 'rows', 'stable');
    p90_summary = groupsummary(citation_track2(uniquerows, ["fcit020", "CPC", "Time"]), ["CPC", "Time"], func);
    p90_summary.p90 = p90_summary.fun1_fcit020; 
    p90_summary = p90_summary(:, ["CPC", "Time", "p90"]);

    citation_track2 = outerjoin(citation_track2, p90_summary, 'Type', 'left', 'MergeKeys', true, 'Keys', ["CPC", "Time"]);
    citation_track2.Major = citation_track2.fcit020 >= citation_track2.p90;


    quality = groupsummary(citation_track2(:, ["CPC", "Firm", "Time", "Major", "count"]), ...
                                ["CPC", "Firm", "Time"], ...
                                {'sum', 'nnz'});
    quality.quality = quality.sum_Major./quality.nnz_count;
   


end