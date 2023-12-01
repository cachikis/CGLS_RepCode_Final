%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: finfric_scatter.m
% Author: Craig A. Chikis
% Date: 10/06/2023
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = finfric_scatter()
m = benchmark_new(); 


beta_vec = linspace(0.015, 0.023, 8)/12; 

mcell = cell(1, length(beta_vec)); 
for ii = 1:length(beta_vec)
    m = robust_multistep_1(); 
    m.alpha_fric = beta_vec(ii);

    m = workhorse_nested_robust(m); 
    m = sim_wrap(m);
    m = CSTAT_sim(m);  
    mcell{ii} = m; 
end


corr_res =cell(size(mcell)); 
xvec = zeros(size(mcell, 2), 2*mcell{1}.s_max+1); 
muvec = zeros(size(xvec, 1), mcell{1}.s_max+1); 
for ii = 1:length(mcell)
    panel_tmp = mcell{ii}.panel_out;

    panel_tmp.rdgrowth_w = panel_tmp.rdgrowth_w - mean(panel_tmp.rdgrowth_w);
    panel_tmp.profitgrowth_w = panel_tmp.profitgrowth_w - mean(panel_tmp.profitgrowth_w); 


    corr_res{ii} = groupsummary(panel_tmp(:, ["Friction_t", "rdgrowth_w", "profitgrowth_w"]), ["Friction_t"], @(x,y) corr_coef_helper(x,y), ...
                                {["rdgrowth_w", "rdgrowth_w"], ["profitgrowth_w", "profitgrowth_w"]});


    xvec(ii, :) = mcell{ii}.investment; 
    muvec(ii, :) = mcell{ii}.firm_distributions; 

end

sz = cellfun(@(x) x.GroupCount(x.Friction_t == 0), corr_res) / 100;
vec_constrained = zeros(2, size(corr_res, 2)); 
for ii = 1:length(corr_res)
    if length(corr_res{ii}.fun1_rdgrowth_w_profitgrowth_w(corr_res{ii}.Friction_t == 1)) > 0
        vec_constrained(1, ii) = corr_res{ii}.fun1_rdgrowth_w_profitgrowth_w(corr_res{ii}.Friction_t == 1);
        vec_constrained(2, ii) = corr_res{ii}.GroupCount(corr_res{ii}.Friction_t == 1) / 100; 
    else
        vec_constrained(1, ii) = nan; 
        vec_constrained(2, ii) = 1e-4 / 100;  
    end
end


close all 
figure; 


set(gcf, 'PaperUnits', 'inches');
x_width=4;
y_width=3;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


scatter(beta_vec(1:(end))*12, cellfun(@(x) x.fun1_rdgrowth_w_profitgrowth_w(x.Friction_t == 0), ...
        corr_res(1:(end))), ...
        sz, "blue", "filled");
hold on 
scatter(beta_vec(1:(end))*12, vec_constrained(1, 1:(end)), ...
        vec_constrained(2, 1:(end)), "red", "filled");
ylim([-0.1,1])
xlim([0.014, 0.025])
yticks([0,0.5,1])
xlabel('Pledgeable share, $\beta$', 'Interpreter', 'latex')
ylabel('Correlation', 'Interpreter', 'latex')
l1 = legend('Unconstrained', 'Constrained', 'Location', 'east', 'Interpreter', 'latex');
set(l1, 'box', 'off')
box off; 

saveas(gcf, 'Output/Figures_Paper/fric_fact.png')
saveas(gcf, 'Output/Figures_Paper/fric_fact.eps', 'epsc')

end