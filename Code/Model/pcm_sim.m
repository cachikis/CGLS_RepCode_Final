%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: pcm_sim.m
% Author: Craig A. Chikis
% Date: 10/18/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = pcm_sim(m)


    panel = m.panel_save; 
    citation_track2 = m.citation_track2; 
    month_increment = m.month_increment;
    innovation_mat = m.innovation_mat; 
    fcitcount = m.fcitcount;  
    nbucket_patent = m.nbucket_patent; 

    condition_uniform =  ~(panel.xrd > 0 & panel.Sale_t > 0 & panel.Profit_t > 0 & ...
						   ~isnan(panel.rdsales) & ~isnan(panel.profitgrowth) & ~isnan(panel.salegrowth) & ...
						   ~isnan(panel.Profit_t) & ~isnan(panel.Profit_tm1) & ~isnan(panel.Profit_t1));


    panel(condition_uniform, :) = []; 
    panel = bucket_quantile(panel, "Profit_t", nbucket_patent);
    panel = renamevars(panel, ["bucket", "bucket_t1", "bucket_tm1"], ...
                       ["bucket_patent", "bucket_t1_patent", "bucket_tm1_patent"]);
    add_buckets = bucket_quantile(panel, "Profit_t", 2);
    add_buckets = renamevars(add_buckets, ["bucket", "bucket_t1", "bucket_tm1"], ...
                       ["bucket_patent_2", "bucket_t1_patent_2", "bucket_tm1_patent_2"]); 
    panel = outerjoin(panel, add_buckets(:, ["Time", "Firm", "bucket_patent_2", "bucket_t1_patent_2", "bucket_tm1_patent_2"]), ...
                      'Type', 'left', 'MergeKeys', true, 'keys', ["Time", "Firm"]); 
   

    quality = unique(citation_track2(:, ["CPC", "Firm", "Time"]), 'rows');
    quality(quality.Time > max(quality.Time), :) = [];
    quality.quality = nan(size(quality, 1), 1); 
 
    for (ii = 1:size(quality, 1))
            iter = citation_track2(citation_track2.CPC == quality.CPC(ii) & ...
                                   citation_track2.Firm == quality.Firm(ii) & ...
                                   ((citation_track2.Time >= quality.Time(ii)) & (citation_track2.Time < quality.Time(ii) + 1)), :);

            quality.quality(ii) = sum(iter.Major)/sum(~isnan(iter.fcit020)); 
    end
   
    panelmerge = panel;
    panelmerge.Year = panelmerge.Time; 
    qualitysale = innerjoin(panelmerge, quality, 'keys', ["Firm", "Time"]); 
    qualityprof_tm1 = qualitysale;

    qualityprof_tm1.bucket_tm1 = nominal(qualityprof_tm1.bucket_tm1_patent);
    qualityprof_tm1.CPC = nominal(qualityprof_tm1.CPC);
    qualityprof_tm1.Time = nominal(qualityprof_tm1.Time);


    cdf_uncond = zeros(7,2);
    cdf_uncond(:, 1) = [0,3,4,5,10,20,50]';
    cdf_uncond(1,2) = sum(citation_track2.fcit020 <= 0)/sum(~isnan(citation_track2.fcit020));
    cdf_uncond(2,2) = sum(citation_track2.fcit020 <= 3)/sum(~isnan(citation_track2.fcit020));
    cdf_uncond(3,2) = sum(citation_track2.fcit020 <= 4)/sum(~isnan(citation_track2.fcit020));
    cdf_uncond(4,2) = sum(citation_track2.fcit020 <= 5)/sum(~isnan(citation_track2.fcit020));
    cdf_uncond(5,2) = sum(citation_track2.fcit020 <= 10)/sum(~isnan(citation_track2.fcit020));
    cdf_uncond(6,2) = sum(citation_track2.fcit020 <= 20)/sum(~isnan(citation_track2.fcit020));
    cdf_uncond(7,2) = sum(citation_track2.fcit020 <= 50)/sum(~isnan(citation_track2.fcit020));
    cdf_tab = array2table(cdf_uncond, 'VariableNames', {'fcit020', 'Prob'});


    m.cdf_tab = cdf_tab;
    m.qualityprof_tm1 = qualityprof_tm1; 

  

     


end