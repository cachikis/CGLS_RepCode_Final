%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: main_run.m
% Author: Craig A. Chikis
% Date: 08/07/2023
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = main_run()
m = benchmark_new(); 
m = workhorse_nested_robust(m); 
m = sim_wrap(m); 
m = CSTAT_sim(m);
m = innovation_output(m); 

% Numbers for model quantitification
chi = 2.143390;
rk = 5.959426;
rhotilde = rk - 100*((1+m.growth_rate)^12-1) - chi; 


rk_obv = rhotilde + chi + 100*((1+m.growth_rate)^12-1); 
rf_obv = rhotilde - chi + 100*((1+m.growth_rate)^12-1); 


string_main_text = "With these parameter values, we exactly match the targeted return on capital, " + ...
                   num2str(rk_obv, '%0.2f') + "\%, and the targeted risk-free rate, " + num2str(rf_obv, '%0.2f') + "\%. " + ...
                   "In the second step,  we identify the parameters $(\phi, \lambda, B)$ using the simulated method of moments. " + ...  
                   "Comparing moments in the model and the data, we choose the parameters that minimize the criterion.";
writematrix(string_main_text, "Output/LaTeX_Output/string_rk_rf.txt", 'QuoteStrings', false);

run_benchmark(); 
% [status, ~] = system('Rscript Code/Data/outregs_model_data_v4.R');

% Read in targets
profitvoltarg = readtable("Output/Store_Data/profitvol_t_three_dropm1.csv");
rdsaletarg = readtable("Output/Store_Data/rdsales_t_targets_three_dropm1.csv"); 
IOtarg = readtable("Output/Store_Data/kogan_targets.csv"); 

tab1 = "\begin{tabular}{ccc|clcc}  \hline \hline  \multicolumn{2}{c}{Parameter estimates} & & &  " +  ...
       " \multicolumn{3}{c}{Moments used in estimation} \\  \cline{1-2} \cline{5-7}   \multicolumn{1}{l}{Parameter} & " + ...  
       " \multicolumn{1}{c}{Value} & & &   \multicolumn{1}{l}{Description} &  \multicolumn{1}{c}{Model} & \multicolumn{1}{c}{Data} \\ " + ...
       " \hline $\phi$ & " + num2str(m.phi_wt, '%0.3f') + " &  &  & Productivity growth & " + num2str(100*((1+m.growth_rate)^12-1), '%0.2f') + "\%" + ...
       " & 1.03\% \\  $\lambda$ & " + num2str(m.lambda, '%0.3f') + " & & &  Mean markup & " + num2str(100*m.markup, '%0.2f') + "\%" + ...
       " & 19.41\% \\  $B$ & " + num2str(m.B*12, '%0.3f') + " & & & Profit volatility & & \\  & & & & \hspace{0.2in} All firms & " + ...
       num2str(100*m.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + "\%" + " & " + num2str(100*profitvoltarg.sd(end), '%0.2f') + ...
       "\% \\  & & & & \hspace{0.2in} Top profit quintile & " + ...
       num2str(100*m.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + "\%" + " & " + num2str(100*profitvoltarg.sd(1), '%0.2f') + "\% \\ " + ...
       " & & & & R\&D to sales \\ " + ...
       " & & & & \hspace{0.2in} All firms & " + num2str(100*m.std_table.mean_fun2_rdsales_w(end), '%0.2f') + "\%" + ...
       " & " + num2str(100*rdsaletarg.p50(end), '%0.2f') + "\% \\ " + ...
       " & & & & \hspace{0.2in} Top profit quintile & " + num2str(100*m.std_table.mean_fun2_rdsales_w(1), '%0.2f') + "\%" + ...
       " & " + num2str(100*rdsaletarg.p50(1), '%0.2f') + "\% \\ " + ...
       " & & & & Innovation output & & \\ & & & & \hspace{0.2in} Mean & " + num2str(100*m.uncond_IO(1), '%0.2f') + "\%" +  " & " + ...
       num2str(IOtarg.mean(1), '%0.2f') +  "\% \\ " + ...
       " & & & & \hspace{0.2in} 90th percentile & " + num2str(m.uncond_IO(9)*100, '%0.2f') + "\%" + ...
       " & " + num2str(IOtarg.p90(1), '%0.2f') + "\% \\ \hline \end{tabular}"; 

writematrix(tab1, "Output/LaTeX_Output/tab1.txt", 'QuoteStrings', false); 


model_reg = readtable(strcat("Output/Store_Data/panel_out_", num2str(m.numCPC), "_", "benchmark", ".csv"));

footnote_regressions = "\footnote{In regressions using Compustat data, an industry is assigned to each firm using the two-digit level of the " + ...
                       "Standard Industrial Classification system.  The regression sample includes " + ...
                       num2str(floor(size(model_reg, 1)/1000), '%0.0f') + "," + num2str(size(model_reg, 1) - 1000*floor(size(model_reg, 1)/1000), '%0.0f') + ...
                       " observations in Panel A and  58,461 observations in Panel B.  " + ...
                       "Profit growth and R\&D-to-sales regressions are calculated for a consistent sample so that both variables are " + ...
                       "defined and the number of observations is equal in the two regressions reported in each panel.}"; 
writematrix(footnote_regressions, "Output/LaTeX_Output/footnote_regressions.txt", 'QuoteStrings', false); 

close all
figure;

uncond_IO_targets = readmatrix(strcat("Output/Store_Data/kogan_targets.csv"));
rng(2022+08+30, 'Threefry')
draw_Lerner = betarnd(1.36, 8, 1, 1e6);
draw_Markup = 1./(1-draw_Lerner);

if (m.kappa >= 9999)
       lerner_p = 1-m.lambda.^(-min(0:m.s_max,m.LMS_bar_s));
       markup_p = m.lambda.^min(0:m.s_max, m.LMS_bar_s);
else
       lerner_p = 1-1./m.markup_wt_array;
       markup_p = m.markup_wt_array; 
end

mid_point_markup = conv(m.lambda.^(min(0:m.s_max, m.LMS_bar_s)), [.5,.5], 'valid');
mid_point_lerner = conv(1-m.lambda.^(-min(0:m.s_max, m.LMS_bar_s)), [.5,.5], 'valid');
frac = zeros(1,m.s_max+1);
frac_l = zeros(1,m.s_max+1);
for (ii = 2:length(mid_point_markup)-1)
       frac(ii) = sum(draw_Markup >= mid_point_markup(ii-1) & (draw_Markup < mid_point_markup(ii)));
       frac_l(ii) = sum(draw_Lerner >= mid_point_lerner(ii-1) & (draw_Lerner < mid_point_lerner(ii)));
end 
frac_l(1) = sum(draw_Lerner < mid_point_lerner(1));
frac_l(end) = sum(draw_Lerner >= mid_point_lerner(end));

frac(1) = sum(draw_Markup < mid_point_markup(1));
frac(end) = sum(draw_Markup >= mid_point_markup(end));
frac = frac/sum(frac);

frac_l = frac_l/sum(frac_l);




close all
figure;

uncond_IO_targets = readmatrix(strcat("Output/Store_Data/kogan_targets.csv"));
rng(2022+08+30, 'Threefry')
draw_Lerner = betarnd(1.36, 8, 1, 1e6);
draw_Markup = 1./(1-draw_Lerner);

xprc = [1,5,10,25,50,75,90,95,99]; 
for ii = 1:length(xprc)
       [model_perc(ii), ~, ~] =  markup_quant(m.kappa, m.lambda, m.s_max, m.firm_distributions, m.nu_s, xprc(ii)/100);
end
model_perc = 100*(model_perc - 1);

set(gcf, 'PaperUnits', 'inches');
x_width=1.1*6.5;
y_width=1.1*2.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

subplot(1,2,2);
p1 = plot(xprc, uncond_IO_targets(3:end), '-ob', 'LineWidth', 2);
hold on ;
p2 = plot(xprc, m.uncond_IO(3:end)*100, '--or', 'LineWidth', 2);
title({'Innovation output distribution'; ''}, 'Interpreter', 'latex');
xlabel('Percentile', 'Interpreter', 'latex');
ylabel('Innovation output (\%)', 'Interpreter', 'latex');


subplot(1,2,1);
p2 = plot(xprc, model_perc, '--or', 'LineWidth', 2);
hold on ;
p1 =plot(xprc, 100*(prctile(draw_Markup, xprc) - 1), '-ob', 'LineWidth', 2);
title({'Markup distribution', ''}, 'Interpreter', 'latex');
xlabel('Percentile', 'Interpreter', 'latex');
ylabel('Net markup (\%)', 'Interpreter', 'latex');
l3 = legend([p1, p2], 'Data', 'Model', 'Interpreter', 'latex', 'Location', 'northwest');
set(l3, 'box', 'off');

saveas(gcf, strcat("Output/Figures_Paper/fig1_a.eps"), 'epsc');
saveas(gcf, strcat("Output/Figures_Paper/fig1_a.png"));

% Patent modeling
m = significance_construct_sim(m); 
m = pcm_sim(m);
cdf_data = readtable("Output/Store_Data/kpss1.csv");
nber_data = readtable("Output/Store_Data/nber_cit.csv"); 

[~, idx] = unique(cdf_data.wku);
cdf_data = cdf_data(idx, :); 

[~, idx] = unique(nber_data.wku);
nber_data = nber_data(idx, :); 


close all 
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=3.25;
y_width=2.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


set(gcf, 'PaperUnits', 'inches');
x_width=3.25;
y_width=2.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


xprc = [0,5,10,25,50,100];

data_cdf = zeros(size(xprc));
data_cdf_nber = zeros(size(xprc));
for ii = 1:length(xprc)
       data_cdf(ii) = sum(cdf_data.fcit05_adj <= xprc(ii)) / sum(~isnan(cdf_data.fcit05_adj));
       data_cdf_nber(ii) = sum(nber_data.fcit05_adj <= xprc(ii)) / sum(~isnan(nber_data.fcit05_adj)); 
end

p1 = plot(xprc, 100*data_cdf, '-ok', 'LineWidth', 2);
hold on 
p2 = plot(xprc, 100*data_cdf_nber, '--or', 'LineWidth', 2);
xlabel('Number of citations received', 'Interpreter', 'latex');
ylabel('Cumulative probability', 'Interpreter', 'latex')
l1 = legend('KPST', 'NBER', 'Interpreter', 'latex', 'Location', 'southeast');
set(l1, 'box', 'off')

saveas(gcf, strcat("Output/Figures_Paper/citdist_compare", ".png"));
saveas(gcf, strcat("Output/Figures_Paper/citdist_compare", ".eps"), 'epsc');



close all
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=3.25;
y_width=2.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


xprc = [0,5,10,25,50,100];

data_cdf = zeros(size(xprc));
for ii = 1:length(xprc)
       data_cdf(ii) = sum(cdf_data.fcit05_adj <= xprc(ii)) / sum(~isnan(cdf_data.fcit05_adj));
       model_cdf(ii) = sum(m.citation_track2.fcit020 <= xprc(ii)) / sum(~isnan(m.citation_track2.fcit020)); 
end

p1 = plot([0, xprc], [0, data_cdf*100], '-ob', 'LineWidth', 2); 
hold on 
p2 = plot([0, xprc], [0, model_cdf*100], '--or', 'LineWidth', 2);
xlabel('Number of citations received', 'Interpreter', 'latex');
ylabel('Cumulative probability', 'Interpreter', 'latex');
l1 = legend([p1,p2], 'Data', 'Model', 'Interpreter', 'latex', 'Location', 'southeast');
set(l1,'box','off');

saveas(gcf, strcat("Output/Figures_Paper/citdist_cdf", ".png"));
saveas(gcf, strcat("Output/Figures_Paper/citdist_cdf", ".eps"), 'epsc');


cit_table = strcat("\begin{tabular}{ccc|clcc}  \hline \hline \multicolumn{2}{c}{Parameter estimate} & & & ", ... 
                   " \multicolumn{3}{c}{Moments used in estimation} \\ \cline{1-2} \cline{5-7}  ", ...
                    " \multicolumn{1}{l}{Parameter} &  \multicolumn{1}{c}{Value} & & &  ", ...
                    " \multicolumn{1}{l}{Description} &  \multicolumn{1}{c}{Model} & \multicolumn{1}{c}{Data} \\", ...
                    " \hline & & & & \\ $\vartheta_p$ & ", num2str(m.zeta, '%0.3f'), " & ", ...
                    " & & $\Pr\{\text{Cit.} \leq 5\}$ & ", num2str(100*model_cdf(2), '%0.2f'), "\%", " & ", ...
                            num2str(100*data_cdf(2), '%0.2f'), "\% \\ & & & & \\ ", ...
                    " & & & & $\Pr\{\text{Cit.} \leq 10\}$ & ", num2str(100*model_cdf(3), '%0.2f'), "\%", " & ", ...
                            num2str(100*data_cdf(3), '%0.2f'), "\% \\ ", ...
                      "\hline \end{tabular}");
writematrix(cit_table, "Output/LaTeX_Output/citation_table.txt", 'QuoteStrings', false); 



close all 
figure;

m = transmat_sim(m); 
m = FHK_sim(m); 


profvol_targets = readtable(strcat("Output/Store_Data/profitvol_t_three_dropm1", ".csv")); 

data_transmat_main = readmatrix(strcat("Output/Store_Data/transmat_targets_", "three_dropm1", ".csv")); 
data_transmat_top10 = readmatrix(strcat("Output/Store_Data/transmat_targets_bottom90_", "three_dropm1.csv"));
data_transmat_top20 = readmatrix(strcat("Output/Store_Data/transmat_targets_top10top20_", "three_dropm1.csv"));
data_transmat_halves = readmatrix(strcat("Output/Store_Data/transmat_targets_halves_", "three_dropm1", ".csv")); 

bottom60top20 = sum(data_transmat_main(3:5, 1)); 
top10bottom90 = sum(data_transmat_top10(1,2));
top50bottom50 = sum(data_transmat_halves(1,2)); 
top20top20 = sum(data_transmat_top20(1,1)); 

bottom60top20_model = sum(m.transmat_out(3:5, 1)); 
top10bottom90_model = sum(m.transmat_out_top10(1,2));
top50bottom50_model = sum(m.transmat_out_halves(1,2));
top20top20_model = sum(m.transmat_out_top20(1,1));

total_prod = 3.44 + -.41 + .76 + 1.23 + .12; 
ENTRANCE_data = 1.23/total_prod;
WITHIN_data = 3.44/total_prod;
BETWEEN_data = -.41/total_prod;
CROSS_data = .76/total_prod;
EXIT_data = .12/total_prod;
    
FHK_targ = 100*[WITHIN_data*total_prod/(total_prod - ENTRANCE_data*total_prod - EXIT_data*total_prod), ...
                    BETWEEN_data*total_prod/(total_prod - ENTRANCE_data*total_prod - EXIT_data*total_prod), ...
                    CROSS_data*total_prod/(total_prod - ENTRANCE_data*total_prod - EXIT_data*total_prod), ...
                    1e-8,1e-8];


set(gcf, 'PaperUnits', 'inches');
x_width=8;
y_width=3.1;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);

subplot(1,3,1);
X = categorical({'Q1','Q2','Q3','Q4', 'Q5'});
X = reordercats(X,{'Q1','Q2','Q3','Q4', 'Q5'});
b1 = bar(X, 100*[profvol_targets.sd(1:5), m.std_table.mean_fun1_profitgrowth_w(1:5)], 'FaceColor', 'flat');
b1(1).CData = [0,0,1];
b1(2).CData = [1,0,0]; 
title({'Profit volatility'; ''}, 'Interpreter', 'latex');
ax = gca;
ax.XAxis.TickLabelInterpreter= 'latex';
ax.XAxis.FontSize = 8;
xlabel('Size quintile, Q5 is smallest quintile', 'Interpreter', 'latex');
l1 = legend('Data', 'Model', 'Interpreter', 'latex', 'Location', 'northwest');
set(l1, 'box', 'off');

subplot(1,3,2)
X = categorical({'Top 50\% to Bottom 50\%', 'Top 10\% to Bottom 90\%', 'Top 20\% to Top 20\%', 'Bottom 60\% to Top 20\%'});
X = reordercats(X,{'Top 50\% to Bottom 50\%', 'Top 10\% to Bottom 90\%', 'Top 20\% to Top 20\%', 'Bottom 60\% to Top 20\%'});
b1 = bar(X, 100*[top50bottom50, top50bottom50_model;
                           top10bottom90, top10bottom90_model;
                              top20top20, top20top20_model;
                              bottom60top20, bottom60top20_model], ...
              'FaceColor', 'flat');
b1(1).CData = [0,0,1];
b1(2).CData = [1,0,0]; 
ax = gca;
ax.XAxis.TickLabelInterpreter= 'latex';
ax.XAxis.FontSize = 8;

title({'Transition rates'; ''}, 'Interpreter', 'latex');

subplot(1,3,3);
X = categorical({'WITHIN', 'BETWEEN', 'CROSS'});
X = reordercats(X, {'WITHIN', 'BETWEEN', 'CROSS'}); 
b1 = bar(X, [FHK_targ(1), 100*m.WITHIN_out;
                      FHK_targ(2), 100*m.BETWEEN_out;
                       FHK_targ(3), 100*m.CROSS_out], 'FaceColor', 'flat');
b1(1).CData = [0,0,1];
b1(2).CData = [1,0,0]; 
title({'FHK decomposition'; ''}, 'Interpreter', 'latex');
ax = gca;
ax.XAxis.TickLabelInterpreter= 'latex';
ax.XAxis.FontSize = 8;

saveas(gcf, strcat("Output/Figures_Paper/2x2_new", ".eps"), 'epsc');
saveas(gcf, strcat("Output/Figures_Paper/2x2_new", ".png"));


% Transition rates text in main
transition_rates_string = "Figure \ref{fig:nontargeted}, middle panel, shows transition rates. Following \citet{Acemoglu_et_al_AER_2018}, " + ...
                          "we start by focusing on the transition rate from the larger half of firms to the smaller half.  In the data, about " + ...
                          num2str(top50bottom50*100,'%0.2f') + "\% of firms in the top half each year are ranked in the bottom half the following year, " + ...
                          "close to the value in the model, " + num2str(100*top50bottom50_model, '%0.2f') + "\%.  We also take a ``best versus the rest'' " + ...
                          "perspective and examine exit from the top decile of firms.  In the data, only " + num2str(100*top10bottom90, '%0.2f') + ... 
                          "\% of firms in the top decile each year are ranked in the bottom 90\% the following year, also close to the value in the model.  " + ...
                          "In the model and the data, firms from the top 20\% mostly remain in the top quintile the following year, " + ...
                          "while almost no firms from the bottom 60\% reach the top quintile.  Thus, our model matches well salient facts about the " + ...
                          "persistence and attainment of leadership.";
writematrix(transition_rates_string, "Output/LaTeX_Output/transition_rate.txt", 'QuoteStrings', false); 



chi = 2.143390;
rk = 5.959426;
rhotilde = rk - 100*((1+m.growth_rate)^12-1) - chi; 
chi_vec = sort([linspace(0, 0.0029, 50), (1+chi/100)^(1/12)-1]); 
rhotilde_vec = sort([linspace(8.46e-4, 0.004, 50), (1+rhotilde/100)^(1/12)-1]);
rfvec = zeros(1, length(chi_vec));
gvec = zeros(1, length(chi_vec)); 
mcell = cell(1, length(chi_vec));
rfvec_rhotilde = zeros(1, length(chi_vec));
gvec_rhotilde = zeros(1, length(chi_vec)); 
markupvec = zeros(size(gvec));
markupvec_rhotilde = zeros(size(gvec_rhotilde)); 
logtfprdispersion_vec = zeros(size(gvec_rhotilde)); 
logtfprdispersion_vec_rhotilde = zeros(size(gvec_rhotilde)); 
emp_growth_dispersion_vec = zeros(size(gvec_rhotilde)); 
emp_growth_dispersion_vec_rhotilde = zeros(size(gvec_rhotilde));
emp_growth_gross_reallocation_vec = zeros(size(gvec_rhotilde)); 
emp_growth_gross_reallocation_vec_rhotilde = zeros(size(gvec_rhotilde));
hhi_vec = zeros(size(gvec_rhotilde)); 
hhi_vec_rhotilde = zeros(size(gvec_rhotilde));  
profitshare_vec = zeros(size(gvec_rhotilde));
profitshare_vec_rhotilde = zeros(size(gvec_rhotilde)); 
for ii = 1:length(chi_vec)
       miter = m; 
       miter.rho = (1+rhotilde/100)^(1/12)-1 + chi_vec(ii); 
       mcell{ii} = miter;
end
parfor ii = 1:length(mcell)
       miter = workhorse_nested_robust(mcell{ii}); 
       miter = sim_wrap(miter); 
       miter = CSTAT_sim(miter); 
       miter = hhi(miter);
       miter = logtfprdispersion(miter); 
       miter = emp_growth_dispersion(miter); 
       miter = profitshare(miter); 

       gvec(ii) = 100*((1+miter.growth_rate)^12-1); 
       rfvec(ii) = 100*((1 + (1+rhotilde/100)^(1/12)-1 - chi_vec(ii) + miter.growth_rate)^12-1); 
       markupvec(ii) = 100*miter.markup; 
       logtfprdispersion_vec(ii) = miter.logtfprdispersion; 
       emp_growth_dispersion_vec(ii) = miter.emp_growth_dispersion;
       emp_growth_gross_reallocation_vec(ii) = miter.emp_growth_gross_reallocation;
       hhi_vec(ii) = miter.hhi; 
       profitshare_vec(ii) = miter.profitshare; 
end

for ii = 1:length(chi_vec)
       miter = m; 
       miter.rho = rhotilde_vec(ii) + (1+chi/100)^(1/12)-1; 
       mcell{ii} = miter;
end
parfor ii = 1:length(mcell)
       miter = workhorse_nested_robust(mcell{ii}); 
       miter = sim_wrap(miter); 
       miter = CSTAT_sim(miter); 
       miter = hhi(miter);
       miter = logtfprdispersion(miter); 
       miter = emp_growth_dispersion(miter); 
       miter = profitshare(miter); 

       gvec_rhotilde(ii) = 100*((1+miter.growth_rate)^12-1); 
       rfvec_rhotilde(ii) = 100*((1 + rhotilde_vec(ii) - (1+chi/100)^(1/12)+1 + miter.growth_rate)^12-1); 
       markupvec_rhotilde(ii) = 100*miter.markup;
       gvec(ii) = 100*((1+miter.growth_rate)^12-1); 
       rfvec(ii) = 100*((1 + (1+rhotilde/100)^(1/12)-1 - chi_vec(ii) + miter.growth_rate)^12-1); 
       markupvec(ii) = 100*miter.markup; 
       logtfprdispersion_vec_rhotilde(ii) = miter.logtfprdispersion; 
       emp_growth_dispersion_vec_rhotilde(ii) = miter.emp_growth_dispersion;
       emp_growth_gross_reallocation_vec_rhotilde(ii) = miter.emp_growth_gross_reallocation;
       hhi_vec_rhotilde(ii) = miter.hhi; 
       profitshare_vec_rhotilde(ii) = miter.profitshare; 
end

close all 
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=7.75;
y_width=2.75;
set(gcf, 'PaperPosition', [0 0 x_width y_width]);

subplot(1,3,1)
plot(rfvec, 100*emp_growth_dispersion_vec, '-b', 'LineWidth', 2);
hold on 
plot(rfvec_rhotilde, 100*emp_growth_dispersion_vec_rhotilde, '--r', 'LineWidth', 2);
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex')
ylabel('\%', 'Interpreter', 'latex')
title({'Firm growth dispersion'; ''}, 'Interpreter', 'latex')

subplot(1,3,2)
plot(rfvec, 100*emp_growth_gross_reallocation_vec, '-b', 'LineWidth', 2);
hold on 
plot(rfvec_rhotilde, 100*emp_growth_gross_reallocation_vec_rhotilde, '--r', 'LineWidth', 2);
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex')
ylabel('\%', 'Interpreter', 'latex')
title({'Gross reallocation'; ''}, 'Interpreter', 'latex')
l1 = legend('Varying $\chi$', 'Varying $\rho$', 'Interpreter', 'latex', 'Location', 'north');
ylim([15, 21])
set(l1, 'box', 'off')

subplot(1,3,3)
plot(rfvec, 100*profitshare_vec, '-b', 'LineWidth', 2);
hold on 
plot(rfvec_rhotilde, 100*profitshare_vec_rhotilde, '--r', 'LineWidth', 2);
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex')
ylabel('\%', 'Interpreter', 'latex')
title({'Profit share'; ''}, 'Interpreter', 'latex')


saveas(gcf, "Output/Figures_Paper/additional_moments.png")
saveas(gcf, "Output/Figures_Paper/additional_moments.eps", 'epsc')






close all 
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=6.5;
y_width=2.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(1,2,1)
p1 = plot(linspace(0,max(rfvec),5), spline(rfvec, gvec, linspace(0, max(rfvec), 5)), '-b', 'LineWidth', 2);
hold on 
p2 = plot(rfvec_rhotilde, gvec_rhotilde, '--r', 'LineWidth', 2);
xlabel('Risk-free rate', 'Interpreter', 'latex');
title({'Growth'; ''}, 'Interpreter', 'latex')
ylabel('\%', 'Interpreter', 'latex')
ylim([0.9, 1.25]);
l1= legend([p1,p2], 'Varying $\chi$', 'Varying $\rho$', 'Interpreter', 'latex', 'Location', 'northeast');
set(l1, 'box', 'off');



subplot(1,2,2)
p1 = plot(linspace(0,max(rfvec),5), spline(rfvec, markupvec, linspace(0,max(rfvec),5)), '-b', 'LineWidth', 2);
hold on 
p2 = plot(rfvec_rhotilde, markupvec_rhotilde, '--r', 'LineWidth', 2);
xlabel('Risk-free rate', 'Interpreter', 'latex');
ylabel('\%', 'Interpreter', 'latex')
title({'Mean net markup'; ''}, 'Interpreter', 'latex');

saveas(gcf, strcat("Output/Figures_Paper/blackboard_new", ".eps"), 'epsc');
saveas(gcf, strcat("Output/Figures_Paper/blackboard_new", ".png"));


% BGP figs
chivec = sort([(1+chi / 100)^(1/12)-1, (1+(5)/100)^(1/12)-1]);
miter = m; 
mcell = cell(1, length(chivec)); 
for ii = 1:length(chivec)
       miter.rho = (1+rhotilde/100)^(1/12)-1 + chivec(ii);
       mcell{ii} = workhorse_nested_robust(miter);
end

close all 
figure;

set(gcf, 'PaperUnits', 'inches');
x_width=6.5;
y_width=2.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(1,2,1)

p1 = plot(0:m.s_max, 12*mcell{1}.investment((m.s_max+1):end), '--b', 'LineWidth', 2);
hold on 
p2 = plot(0:m.s_max, 12*flip(mcell{1}.investment(1:(m.s_max+1))), ':b', 'LineWidth', 2);

p3 = plot(0:m.s_max, 12*mcell{2}.investment((m.s_max+1):end), '--r', 'LineWidth', 2);
hold on 
p4 = plot(0:m.s_max, 12*flip(mcell{2}.investment(1:(m.s_max+1))), ':r', 'LineWidth', 2);

xlabel('Technology gap, $s$','Interpreter', 'latex');
title({'Cross-section of firm innovation rates, $x_{|\sigma|}$'; ''}, 'Interpreter', 'latex');
xlim([0, 30]);
ylim([0.1, 0.6])

l1 = legend(strcat("Leader, $\chi = $", num2str(round(chivec(1)*1200), '%0.0f'), "\%"), ...
            strcat("Laggard, $\chi = $", num2str(round(chivec(1)*1200), '%0.0f'), "\%"), ...
            strcat("Leader, $\chi = $", num2str(round(chivec(2)*1200), '%0.0f'), "\%"), ...
            strcat("Laggard, $\chi = $", num2str(round(chivec(2)*1200), '%0.0f'), "\%"), ...
            'Interpreter', 'latex') ;
set(l1, 'box', 'off');


subplot(1,2,2);

p1 = plot(0:m.s_max, mcell{1}.firm_distributions, '-b', 'LineWidth', 2);
hold on 
p3 = plot(0:m.s_max, mcell{2}.firm_distributions, '-r', 'LineWidth', 2);

xlabel('Technology gap, $s$','Interpreter', 'latex');
title({'Distribution of technology gaps, $\mu_s$'; ''}, 'Interpreter', 'latex');
xlim([0, 30]);
l1 = legend(strcat("Leader, $\chi = $", num2str(round(chivec(1)*1200), '%0.0f'), "\%"), ...
            strcat("Laggard, $\chi = $", num2str(round(chivec(2)*1200), '%0.0f'), "\%"), ...
            'Interpreter', 'latex') ;
set(l1, 'box', 'off');



saveas(gcf, strcat("Output/Figures_Paper/bgp_obj", ".eps"), 'epsc');
saveas(gcf, strcat("Output/Figures_Paper/bgp_obj", ".png"));


m = benchmark_new(); 
m = workhorse_nested_robust(m);
m2 = benchmark_new();
m2.rho = m.rho - ((1+chi/100)^(1/12)-1) + ((1+chi/100)^(1/12)-1) + ((1.02)^(1/12)-1); 
m2 = workhorse_nested_robust(m2);


s1 = "With a rise in the excess return parameter $\chi$ from " + num2str(chi, '%0.2f') + ...
       "\% to " + num2str(chi+2,'%0.2f') + "\%, productivity growth falls from " + num2str(100*((1+m.growth_rate)^12-1), '%0.2f') + ...
       "\% to " + num2str(100*((1+m2.growth_rate)^12-1), '%0.2f') + ...
       "\%, reflecting an intensive margin effect (lower innovation rates conditional on a firm's " + ...
       "position) and an extensive margin effect (a smaller share of high-R\&D, competitive industries)"; 
writematrix(s1, "Output/LaTeX_Output/num1.txt", 'QuoteStrings', false);

end