%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: entry_new_run.m
% Author: Craig A. Chikis
% Date: 08/07/2023
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = entry_new_run() 


m = entry_new(); 
m = workhorse_nested_robust(m); 
m = sim_wrap(m); 
m = CSTAT_sim(m);
m = innovation_output(m); 
m = FHK_sim(m); 
m = emp_share_sim(m); 
[status, ~] = system('Rscript Code/Data/outregs_model_data_v4.R');
[status, ~] = system('Rscript Code/Data/bds.R');


entry_rate = 1200*0.5*sum(m.firm_distributions .* m.investment_entrant); 

m2 = entry_new(); 
m2.rho = m2.rho + (1.02)^(1/12)-1; 
m2 = workhorse_nested_robust(m2); 

entry_rate_2 = 1200*0.5*sum(m2.firm_distributions .* m2.investment_entrant);



% Read in targets
profitvoltarg = readtable("Output/Store_Data/profitvol_t_three_dropm1.csv");
rdsaletarg = readtable("Output/Store_Data/rdsales_t_targets_three_dropm1.csv"); 
IOtarg = readtable("Output/Store_Data/kogan_targets.csv"); 
entry_rate_data = readtable("Output/Store_Data/entry_rate_summ.csv");
bds1 = readtable("Output/Store_Data/bds1.csv");
bds3 = readtable("Output/Store_Data/bds3.csv");


string_main = "We jointly estimate the parameters $(\phi,\lambda,B,\phi_E,B_E,l_E),$ setting the catch-up parameter equal for entrants and laggards " + ...
               "($\phi=\phi_E$).  To identify the entry parameters, we include in the estimation the entry rate and the employment shares of firms " + ...
               "less than 2 and 4 years old. Table \ref{tab:entry_param} shows that the model fits well the aggregate growth rate, average markup, and new " + ...
               "entry moments.  Internet Appendix \ref{app:addl_results} shows that the quality of fit is high as well for the remaining moments.  " + ...
               "Innovating entrants do not leapfrog ($l_E=0$) and catch up to the leader with a " + num2str(100*m.phi_wt_tilde, '%0.0f') + ...
               "\% chance.  Figure \ref{fig:entry_FHK} shows that the model fits well the growth decomposition.  In the estimated model, variation in " + ...
               "excess returns implies low interest rates are associated with low entry rates, weak productivity growth, and higher markups " + ...
               "(Figure \ref{fig:entry_g_rho}).  An increase in the excess return parameter $\chi$ from 2.15\% to 4.15\% implies a decline in the " + ...
               "entry rate from " + num2str(entry_rate, '%0.1f') + "\% to " + num2str(entry_rate_2, '%0.1f') +  ...
               "\%, consistent with the observed trend \citep{AkcigitAtesAEJMacro}."; 
writematrix(string_main, "Output/LaTeX_Output/entry_rate_decline.txt", 'QuoteStrings', false); 



s1 = "\begin{tabular}{lcc}  \hline\hline Moments & \shortstack{Model \\ with entry} & Data  \\ " + ...
     " \hline    Productivity growth & " + num2str(100*((1+m.growth_rate)^12-1), '%0.2f') +  "\% & " + "1.03\% \\ " + ...
    " Mean markup & " + num2str(m.markup*100, '%0.2f') + "\% & 19.41\% \\ " + ...
    " Entry rate & " + num2str(entry_rate, '%0.2f') + "\% & " + num2str(entry_rate_data.median(1), '%0.2f') + "\% \\ " + ...
    " Employment share, $< 2$ years & " + num2str(m.share1*100, '%0.2f') + "\% & " + num2str(bds1.median(1)*100, '%0.2f') + "\%  \\ " + ...
    " Employment share, $< 4$ years & " + num2str(m.share3*100, '%0.2f') + "\% & " + num2str(bds3.median(1)*100, '%0.2f') + "\% \\  " + ...
    " \hline \hline Parameters & &  \\  " + ...
    " \hline  $\phi = \phi_E$ & " + num2str(m.phi_wt, '%0.3f') + ...
    " & \\  $\lambda$ & " + num2str(m.lambda, '%0.3f') + ...
    " & \\  $B$ & " + num2str(m.B*12, '%0.3f') + ...
    " & \\  $B_E$ & " + num2str(m.B_entrant*12, '%0.3f') + ...
    " & \\  $l_E$ & " + num2str(m.l_tilde, '%0.0f') + " & \\  \hline  \end{tabular} "; 

s2 = "\begin{tabular}{lcc} \hline\hline Moments & Model with entry & Data \\ " + ...
     " \hline Profit volatility \\ " + ...
     " \hspace{0.2in} All firms & " + num2str(100*m.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + "\% & " + num2str(profitvoltarg.sd(end)*100, '%0.2f') + ...
     "\% \\" + ...
     " \hspace{0.2in} Top profit quintile & " + num2str(100*m.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') +  "\% & " + ...
     num2str(100*profitvoltarg.sd(1), '%0.2f') + "\% \\" + ...
     " R\&D to sales  \\ " + ...
     " \hspace{0.2in} All firms & " + num2str(100*m.std_table.mean_fun2_rdsales_w(end), '%0.2f') +  "\% & " + ...
     num2str(rdsaletarg.p50(end)*100, '%0.2f') + "\% \\" + ...
     " \hspace{0.2in} Top profit quintile & " + num2str(100*m.std_table.mean_fun2_rdsales_w(1), '%0.2f') +  "\% & " + ...
     num2str(100*rdsaletarg.p50(1), '%0.2f') + "\% \\" + ...
     " Innovation output \\" + ...
     " \hspace{0.2in} Mean & " + num2str(100*m.uncond_IO(1), '%0.2f') + "\% & " + num2str(IOtarg.mean(1), '%0.2f') + "\% \\ " + ...
     " \hspace{0.2in} 90\textsuperscript{th} percentile & " + num2str(100*m.uncond_IO(end-2), '%0.2f') +  "\% & " + ...
     num2str(IOtarg.p90(1), '%0.2f') + "\% \\ " + ...
     " \hline \end{tabular}"; 


writematrix(s1, "Output/LaTeX_Output/entry_table.txt", 'QuoteStrings', false);
writematrix(s2, "Output/LaTeX_Output/entry_table_IA.txt", 'QuoteStrings', false);


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
FHK_targ_ent = 100*[WITHIN_data, BETWEEN_data, CROSS_data, ENTRANCE_data, EXIT_data];

ent_string = categorical(["WITHIN", "BETWEEN", "CROSS", "ENTRY", "EXIT"]); 
ent_string = reordercats(ent_string, ["WITHIN", "BETWEEN", "CROSS", "ENTRY", "EXIT"]);
ent_model = [m.WITHIN_out, m.BETWEEN_out, m.CROSS_out, m.ENTRY_out, m.EXIT_out]; 
    


close all 
figure;

set(gcf, 'PaperUnits', 'inches');
x_width=3.25;
y_width=3;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


b1 = bar(ent_string, [FHK_targ_ent', 100*ent_model'], 'FaceColor', 'flat');
b1(1).CData = [0,0,1];
b1(2).CData = [1,0,0]; 
ylabel('Contribution to TFP (\%)', 'Interpreter', 'latex');
ylim([-20, 100]);
l1 = legend("Data", "Model", 'Interpreter', 'latex');
set(l1, 'box', 'off');
ax = gca;
ax.XAxis.TickLabelInterpreter= 'latex';
ax.XAxis.FontSize = 8;


saveas(gcf, "Output/Figures_Paper/ent_decomp.eps", 'epsc');
saveas(gcf, "Output/Figures_Paper/ent_decomp.png");

chi = 2.143390;
rk = 5.959426;
mtmp = benchmark_new();
mtmp = workhorse_nested_robust(mtmp); 
rhotilde = rk - chi  - 100*((1+mtmp.growth_rate)^12-1); 

chivec = repmat(sort([linspace(0, 0.003, 25), (1+chi/100)^(1/12)-1]), [1,1]);
m_cell = cell(1, length(chivec));
count = 0;
for entry_type = ["Undirected"]
    for ii = 1:(length(chivec)/1) 
        count = count + 1; 
        mtmp = m;
        mtmp.rho = chivec(ii) + (1+rhotilde/100)^(1/12)-1;
        mtmp.entry_type = entry_type; 
        m_cell{count} = mtmp;
    end
end
gvec = zeros(1, length(m_cell)); 
rfvec = zeros(size(gvec)); 
markupvec = zeros(size(gvec)); 
entryratevec = zeros(size(gvec));
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec(ii) = 100*mtmp.markup; 
    rfvec(ii) = (1+rhotilde/100)^(1/12)-1 - chivec(ii) + mtmp.growth_rate; 
    rfvec(ii) = 100*((1+rfvec(ii))^12-1);
    entryratevec(ii) = 1200 * 0.5 * sum(mtmp.firm_distributions .* mtmp.investment_entrant); 

end



chivec = repmat(sort([linspace(8.46e-4, 0.004, 25), (1+chi/100)^(1/12)-1]), [1,1]);
m_cell = cell(1, length(chivec));
count = 0;
for entry_type = ["Undirected"]
    for ii = 1:(length(chivec)/1) 
        count = count + 1; 
        mtmp = m;
        mtmp.rho = (1+chi/100)^(1/12)-1 +  chivec(ii); 
        mtmp.entry_type = entry_type; 
        m_cell{count} = mtmp;
    end
end
gvec_2 = zeros(1, length(m_cell)); 
rfvec_2 = zeros(size(gvec)); 
markupvec_2 = zeros(size(gvec)); 
entryratevec_2 = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec_2(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec_2(ii) = 100*mtmp.markup; 
    rfvec_2(ii) = -(1+chi/100)^(1/12)+1 + chivec(ii) + mtmp.growth_rate; 
    rfvec_2(ii) = 100*((1+rfvec_2(ii))^12-1);
    entryratevec_2(ii) = 1200 * 0.5 * sum(mtmp.firm_distributions .* mtmp.investment_entrant); 

end




close all
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=8.2;
y_width=2.75;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(1,3,2)
p1 = plot(rfvec(1:(length(chivec)/1)), gvec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), gvec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
hold on 
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('\%', 'Interpreter', 'latex')
title({'Growth'; ''}, 'Interpreter', 'latex');
l1 = legend("Varying $\chi$", ...
            "Varying $\rho$", ...
            'Interpreter', 'latex', 'Location', 'north');
l1.FontSize = l1.FontSize*0.8;
set(l1, 'box', 'off');

subplot(1,3,3)
p1 = plot(rfvec(1:(length(chivec)/1)), markupvec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), markupvec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
hold on 
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
title({'Mean net markup'; ''}, 'Interpreter', 'latex');
ylabel('\%', 'Interpreter','latex')

subplot(1,3,1)
p1 = plot(rfvec(1:(length(chivec)/1)), entryratevec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), entryratevec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
hold on 
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
title({'Entry rate'; ''}, 'Interpreter', 'latex');
ylabel('\%', 'Interpreter', 'latex')

saveas(gcf, "Output/Figures_Paper/ent_g_markup.eps", 'epsc');
saveas(gcf, "Output/Figures_Paper/ent_g_markup.png");




end