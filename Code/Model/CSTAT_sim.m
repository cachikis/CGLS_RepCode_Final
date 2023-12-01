%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: CSTAT_sim.m
% Author: Craig A. Chikis
% Date: 08/27/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = CSTAT_sim(m)

    panel = m.panel;
    winsor_vec = m.winsor_vec;
    
    uniqueyears = unique(panel.Time);
    panel.profitgrowth_w = nan(size(panel, 1), 1);
    panel.rdsales_w = nan(size(panel, 1), 1);
    panel.empgrowth_w = nan(size(panel, 1), 1); 
    panel.rdgrowth_w = nan(size(panel, 1), 1); 
    panel.salegrowth_w= nan(size(panel, 1) ,1); 
    panel.salegrowth_w_halti= nan(size(panel, 1) ,1); 
    panel.rdgrowth_w_halti= nan(size(panel, 1) ,1); 
    panel.profitgrowth_w_halti= nan(size(panel, 1) ,1); 
    for (ii = 1:length(uniqueyears))
        iter = panel(panel.Time == uniqueyears(ii), :);

        profitgrowthiter = winsor(iter.profitgrowth, winsor_vec);       
        rdsalesiter = winsor(iter.rdsales, winsor_vec); 
        if ~m.entrance 
            rdgrowthiter = winsor(iter.rdgrowth, winsor_vec);
        end
        salegrowthiter = winsor(iter.salegrowth, winsor_vec);
        if ~m.entrance
            rdgrowthiter_halti = winsor(iter.rdgrowth_halti, winsor_vec);
            salegrowthiter_halti = winsor(iter.salegrowth_halti, winsor_vec);
            profitgrowthiter_halti = winsor(iter.profitgrowth_halti, winsor_vec); 
        end
        if ~m.entrance
            empgrowthiter = winsor(iter.empgrowth, winsor_vec);
        else
            empgrowthiter = nan(size(iter.profitgrowth)); 
        end

        panel.profitgrowth_w(panel.Time == uniqueyears(ii)) = profitgrowthiter;
        panel.rdsales_w(panel.Time == uniqueyears(ii)) = rdsalesiter; 
        panel.empgrowth_w(panel.Time == uniqueyears(ii)) = empgrowthiter;
        if ~m.entrance
            panel.rdgrowth_w(panel.Time == uniqueyears(ii)) = rdgrowthiter; 
        end
        panel.salegrowth_w(panel.Time == uniqueyears(ii)) = salegrowthiter;
        if ~m.entrance
            panel.rdgrowth_w_halti(panel.Time == uniqueyears(ii)) = rdgrowthiter_halti;
            panel.salegrowth_w_halti(panel.Time == uniqueyears(ii)) = salegrowthiter_halti;
            panel.profitgrowth_w_halti(panel.Time == uniqueyears(ii)) = profitgrowthiter_halti;     
        end
      
    end


    std_table = groupsummary(panel(:, ["Time", "bucket", "profitgrowth_w",  "rdsales_w"]), ...
                             ["Time", "bucket"], ...
                             {@(x) nanstd(x), @(x) nanmedian(x)});
    std_table(:, ["Time", "GroupCount"]) = []; 
    std_table = groupsummary(std_table, "bucket", 'mean');
    std_table = std_table(:, ["bucket", "GroupCount", "mean_fun1_profitgrowth_w", "mean_fun2_rdsales_w"]);

    std_table_uncond = groupsummary(panel(:, ["Time", "profitgrowth_w", "rdsales_w"]), ...
                                    ["Time"], ...
                                    {@(x) nanstd(x), @(x) nanmedian(x)});                                                                    
    std_table(end+1, :) = array2table([0, nan, nanmean(std_table_uncond.fun1_profitgrowth_w), nanmean(std_table_uncond.fun2_rdsales_w)]); 

    
    panel_out = panel; 
       

    m.std_table = std_table;
    m.panel_out = panel_out;
 


end