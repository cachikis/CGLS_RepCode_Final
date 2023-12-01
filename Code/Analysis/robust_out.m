%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: robust_out.m
% Author: Craig A. Chikis
% Date: 08/24/2023
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = robust_out() 

mvec = {@robust_noposRD, @robust_nothreeqtr, @robust_multistep_1, @robust_gamma, @robust_varphi, ...
        @robust_elaslabor, @robust_highmarkup, @robust_eta, @robust_finfric, @robust_kappa12, ...
        @robust_onlyincremental}; 


m_res = cell(size(mvec));
for ii = 1:length(mvec)

    m = mvec{ii}(); 
    m = workhorse_nested_robust(m);
    m = sim_wrap(m);  
    m = CSTAT_sim(m);
    m = innovation_output(m); 

    if ii == 3

        
        % Regression 
        m.panel_out.bucket_tm1 = categorical(m.panel_out.bucket_tm1); 
        m.panel_out.Time = categorical(m.panel_out.Time); 
        
        model = fitlm(m.panel_out, 'rdsales_w ~ bucket_tm1 + Time');

        m.loadings = ...
            model.Coefficients.Estimate(ismember(model.Coefficients.Row, {'bucket_tm1_2', 'bucket_tm1_3', 'bucket_tm1_4', 'bucket_tm1_5'}))';

    end

    if ismember(ii, [7, 11])
        [m.p50_markup, ~, ~] = markup_quant(m.kappa, m.lambda, m.s_max, ...
                                            m.firm_distributions, m.nu_s, 0.5); 
        [m.p75_markup, ~, ~] = markup_quant(m.kappa, m.lambda, m.s_max, ...
                                            m.firm_distributions, m.nu_s, 0.75);
        [m.p90_markup, ~, ~] = markup_quant(m.kappa, m.lambda, m.s_max, ...
                                            m.firm_distributions, m.nu_s, 0.9);  
    end

    m_res{ii} = m; 

end

s1 = "\begin{tabular}{lcccccc} \hline \hline  & \multicolumn{2}{c}{No pos. R\&D} & & " + ...
    " \multicolumn{2}{c}{No 3 qtr. profit} \\   \cline{2-3} \cline{5-6}  Moments & Model & Data & " + ...
    " & Model & Data \\   \hline  " + ...
    " Productivity growth  & " + num2str(100*((1+m_res{1}.growth_rate)^12-1), '%0.2f') + "\% & 1.03\%  & & " + ...
    num2str(100*((1+m_res{2}.growth_rate)^12-1), '%0.2f') + "\% & 1.03\% \\ " + ...
    " Mean markup  & " + num2str(100*m_res{1}.markup, '%0.2f') + "\% & 19.41\% & & " + ...
    num2str(100*m_res{2}.markup, '%0.2f') + "\% & 19.41\% \\ " + ...
    " Profit volatility \\ " + ...
    " \hspace{0.2in} All firms  & " + num2str(100*m_res{1}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + ... 
    "\% & 39.98\%  & & " + num2str(100*m_res{2}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + "\%" +  ...
    " & 40.20\% \\  " + ...
    " \hspace{0.2in} Top profit quintile & " + num2str(100*m_res{1}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + ...
    "\% & 25.52\% & & " + num2str(100*m_res{2}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + ...
    "\% & 19.02\% \\  R\&D to sales \\  " + ...
    " \hspace{0.2in} All firms  & " + num2str(100*m_res{1}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + ...
    "\% & 2.67\%  & & " + num2str(100*m_res{2}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + ...
    "\% & 4.29\% \\ " + ...
    " \hspace{0.2in} Top profit quintile  & " + num2str(100*m_res{1}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + ...
    "\% & 1.78\%  & & " + num2str(100*m_res{2}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + ...
    "\% & 2.67\% \\  Innovation output \\ " + ...
    " \hspace{0.2in} Mean  & " + num2str(100*m_res{1}.uncond_IO(1), '%0.2f') + "\% & 3.43\%  & & " + ...
    num2str(100*m_res{2}.uncond_IO(1), '%0.2f') + "\% & 6.78\% \\  " + ...
    " \hspace{0.2in} 90th percentile  & " + num2str(100*m_res{1}.uncond_IO(end-2), '%0.2f') + "\% & 8.07\% & & " + ...
    num2str(100*m_res{2}.uncond_IO(end-2), '%0.2f') + "\% & 19.54\% \\  " + ...
    " \hline \hline Parameters \\  " + ...
    " \hline $\phi$ & " + num2str(m_res{1}.phi_wt, '%0.3f') + " & & & " + num2str(m_res{2}.phi_wt, '%0.3f') + ...
    " \\  $\lambda$ & " + num2str(m_res{1}.lambda, '%0.3f') + " & & & " + num2str(m_res{2}.lambda, '%0.3f') + ...
    " \\  $B$ & " + num2str(m_res{1}.B*12, '%0.3f') + " & & & " + num2str(m_res{2}.B*12, '%0.3f') + ...
    " \\ \hline \end{tabular}";
writematrix(s1, "Output/LaTeX_Output/moment_alt_table.txt", 'QuoteStrings', false);




s2 = "\begin{tabular}{lcc} \hline \hline  Moments  & Model & Data  \\  " + ...
     " \hline Productivity growth  & " + num2str(100*((1+m_res{3}.growth_rate)^12-1), '%0.2f') + "\% & 1.03\% \\ " + ...
     " Net markup & " + num2str(100*(m_res{3}.markup_wt - 1), '%0.2f') + "\% & 19.41\% \\ Profit volatility \\ " + ...
     " \hspace{0.2in} All firms & " + num2str(100*m_res{3}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + ...
     "\% & 35.72\% \\  " + ...
     " \hspace{0.2in} Top profit quintile & " + num2str(100*m_res{3}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + ...
     "\% & 19.48\% \\  R\&D to sales \\  " + ...
     " \hspace{0.2in} All firms & " + num2str(100*m_res{3}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + ...
     "\% & 4.37\% \\ " + ...
     " \hspace{0.2in} Top profit quintile & " + num2str(100*m_res{3}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + ...
     "\% & 2.64\% \\  R\&D to sales, effect relative to top profit quintile \\  " + ...
     " \hspace{0.2in} Second quintile & " + num2str(100*m_res{3}.loadings(1), '%0.2f') + " & 0.70 \\ " + ... 
     " \hspace{0.2in} Third quintile & " + num2str(100*m_res{3}.loadings(2), '%0.2f') + " & 2.03 \\  " + ...
     " \hspace{0.2in} Fourth quintile & " + num2str(100*m_res{3}.loadings(3), '%0.2f') + " & 3.20 \\ " + ... 
     " \hspace{0.2in} Smallest quintile & " + num2str(100*m_res{3}.loadings(4), '%0.2f') + " & 9.27 \\ " + ...
    " \hline \hline Parameters \\  \hline $\phi_M$ & " + num2str(m_res{3}.phi_wt, '%0.3f') +  " \\  " + ...
    " $\lambda$ & " + num2str(m_res{3}.lambda, '%0.3f') + " \\ " + ... 
    " $B$ & " + num2str(m_res{3}.B*12, '%0.3f') + " \\  " + ...
    " $\kappa$ & " + num2str(m_res{3}.kappa, '%0.3f') + " \\ " + ... 
    " $l_M$ & " + num2str(m_res{3}.l, '%0.0f') + " \\ \hline \end{tabular}";
writematrix(s2, "Output/LaTeX_Output/multistep.txt", 'QuoteStrings', false);



s3 = "\begin{tabular}{lcc} \hline \hline  Moments  & Model & Data  \\ " + ...
     " \hline  Productivity growth & " + num2str(100*((1+m_res{7}.growth_rate)^12-1), '%0.2f') + "\% & 1.03\% " + ...
     " \\  Markup \\  \hspace{0.2in} Mean & " + num2str(100*m_res{7}.markup, '%0.2f') + "\% & 36.50\% \\ " + ...
     " \hspace{0.2in} 50th percentile & " + num2str(100*(m_res{7}.p50_markup-1), '%0.2f') + "\% & 19.70\% \\ " + ...
     " \hspace{0.2in} 75th percentile & " + num2str(100*(m_res{7}.p75_markup-1), '%0.2f') + "\% & 42.20\% \\ " + ...
     " \hspace{0.2in} 90th percentile & " + num2str(100*(m_res{7}.p90_markup-1), '%0.2f') + "\% & 84.80\% \\ " + ...
     " Profit volatility \\  \hspace{0.2in} All firms & " + ...
        num2str(100*m_res{7}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + "\% & 35.72\% \\ " + ...
     "\hspace{0.2in} Top profit quintile & " + num2str(100*m_res{7}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + ...
        "\% & 19.48\% \\  " + ...
     " R\&D to sales \\  \hspace{0.2in} All firms & " + num2str(100*m_res{7}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + ...
     "\% & 4.37\% \\  \hspace{0.2in} Top profit quintile & " + ...
     num2str(100*m_res{7}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + ...
     "\% & 2.64\% \\  " + ...
     " Innovation output \\  \hspace{0.2in} Mean & " + num2str(100*m_res{7}.uncond_IO(1), '%0.2f') + ...
     "\% & 6.78\% \\  " + ...
     " \hspace{0.2in} 90th percentile & " + num2str(100*m_res{7}.uncond_IO(end-2), '%0.2f') + ...
     "\% & 19.54\% \\  \hline \hline Parameters \\  " + ...
     " \hline $\phi$ & " + num2str(m_res{7}.phi_wt, '%0.3f') + ...
     " & \\  $\lambda$ & " + num2str(m_res{7}.lambda, '%0.3f') + ...
     " & \\ $B$ & " + num2str(12*m_res{7}.B, '%0.3f') + " & \\  \hline \end{tabular} "; 
writematrix(s3, "Output/LaTeX_Output/highmarkup.txt", 'QuoteStrings', false);

s4 = "\begin{tabular}{lccccc}  \hline\hline Moments & $\gamma = 0.35$ & $\eta =0.03$ & " + ...
     " $\varphi = 0.50$  & \shortstack{Inelastic\\labor} & Data \\ " + ...
     " \hline Productivity growth & " + num2str(100*((1+m_res{4}.growth_rate)^12-1), '%0.2f') + ...
     "\% & " + num2str(100*((1+m_res{8}.growth_rate)^12-1), '%0.2f') + "\% & " + ...
     num2str(100*((1+m_res{5}.growth_rate)^12-1), '%0.2f') + "\% & " + ...
     num2str(100*((1+m_res{6}.growth_rate)^12-1), '%0.2f') + "\% & 1.03\% \\ " + ... 
     " Mean markup & " + num2str(100*m_res{4}.markup, '%0.2f') + "\% & " + num2str(100*m_res{8}.markup, '%0.2f') + ...
     "\% & " + num2str(100*m_res{5}.markup, '%0.2f') + "\% & " + ...
     num2str(100*m_res{6}.markup, '%0.2f') + "\% & 19.41\% \\ " + ...
     " Profit volatility \\  " + ...
     " \hspace{0.2in} All firms & " + num2str(100*m_res{4}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + ...
     "\% & " + num2str(100*m_res{8}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f')  + ...
     "\% & " + num2str(100*m_res{5}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + ...
     "\% & " + num2str(100*m_res{6}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + ...
     "\% & 35.72\% \\ " + ...
     " \hspace{0.2in} Top profit quintile & " + num2str(100*m_res{4}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + ...
     "\% & " + num2str(100*m_res{8}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + ...
     "\% & " + num2str(100*m_res{5}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + ... 
     "\% & " + num2str(100*m_res{6}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + ...
     "\% & 19.48\% \\  R\&D to sales \\  " + ...
     " \hspace{0.2in} All firms & " + num2str(100*m_res{4}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + ...
     "\% & " + num2str(100*m_res{8}.std_table.mean_fun2_rdsales_w(end), '%0.2f')  + ...
     "\% & " + num2str(100*m_res{5}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + ...
     "\% & " + num2str(100*m_res{6}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + ...
     "\% & 4.37\% \\  " + ...
     " \hspace{0.2in} Top profit quintile & " + num2str(100*m_res{4}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + ...
     "\% & " + num2str(100*m_res{8}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + ...
     "\% & " + num2str(100*m_res{5}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + ...
     "\% & " + num2str(100*m_res{6}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + ...
     "\% & 2.64\% \\  Innovation output \\  " + ...
     " \hspace{0.2in} Mean & " + num2str(100*m_res{4}.uncond_IO(1), '%0.2f') + ...
     "\% & " + num2str(100*m_res{8}.uncond_IO(1), '%0.2f') + ...
     "\% & "  + num2str(100*m_res{5}.uncond_IO(1), '%0.2f') + ...
     "\% & " + num2str(100*m_res{6}.uncond_IO(1), '%0.2f') + ...
     "\% & 6.78\% \\  " + ...
     " \hspace{0.2in} 90th percentile & " + num2str(100*m_res{4}.uncond_IO(end-2), '%0.2f') + ...
     "\% & " + num2str(100*m_res{8}.uncond_IO(end-2), '%0.2f') + ...
     "\% & "  + num2str(100*m_res{5}.uncond_IO(end-2), '%0.2f') + ...
     "\% & " + num2str(100*m_res{6}.uncond_IO(end-2), '%0.2f') + ...
     "\% & 19.54\% \\  " + ...
     " \hline \hline Parameters \\  " + ...
     " \hline $\phi$ & " + num2str(m_res{4}.phi_wt, '%0.3f') + " & " + num2str(m_res{8}.phi_wt, '%0.3f') + ...
     " & " + num2str(m_res{5}.phi_wt, '%0.3f') + " & " + num2str(m_res{6}.phi_wt, '%0.3f') + " \\  " + ...
     " $100 \cdot \ln (\lambda)$ & " + num2str(100*log(m_res{4}.lambda), '%0.3f') + ... 
     " & " + num2str(100*log(m_res{8}.lambda), '%0.3f') + " & " + ...
     num2str(100*log(m_res{5}.lambda), '%0.3f') + " & " + num2str(100*log(m_res{6}.lambda), '%0.3f') + ...
     " \\  $B$ & " + num2str(m_res{4}.B*12, '%0.3f') + " & " + ...
     num2str(m_res{8}.B*12, '%0.3f') + " & " + num2str(m_res{5}.B*12, '%0.3f') + ...
     " & " + num2str(m_res{6}.B*12, '%0.3f') + " \\ \hline \end{tabular}"; 
writematrix(s4, "Output/LaTeX_Output/robustness_main.txt", 'QuoteStrings', false)


s5 = "\begin{tabular}{lcc} \hline \hline Moments & $\beta = 0.1$ & Data \\ " + ...
     " \hline Productivity growth  & " + num2str(100*((1+m_res{9}.growth_rate)^12-1), '%0.2f') + "\% & 1.03\% \\ " + ...
     " Mean markup & " + num2str(100*m_res{9}.markup, '%0.2f') + "\% & 19.41\% \\ " + ...
     " Profit volatility \\  " + ...
     " \hspace{0.2in} All firms & " + num2str(100*m_res{9}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + "\% & 35.72\% \\ " + ...
     " \hspace{0.2in} Top profit quintile & " + num2str(100*m_res{9}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + "\% & 19.48\% \\ " + ...
     " R\&D to sales \\ " + ...
     " \hspace{0.2in} All firms & " + num2str(100*m_res{9}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + "\% & 4.37\% \\ " + ...
     " \hspace{0.2in} Top profit quintile & " + num2str(100*m_res{9}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + "\% & 2.64\% \\ " + ...
     " Innovation output \\  " + ...
     " \hspace{0.2in} Mean & " + num2str(100*m_res{9}.uncond_IO(1), '%0.2f') + "\% & 6.78\% \\ " + ...
     " \hspace{0.2in} 90th percentile & " + num2str(100*m_res{9}.uncond_IO(end-2), '%0.2f') + "\% & 19.54\% \\ " + ...
     " \hline \hline Parameters \\  \hline $\phi$ & " + num2str(m_res{9}.phi_wt, '%0.3f') + ...
     " \\  $\lambda$ & " + num2str(m_res{9}.lambda, '%0.3f') + ...
     " \\  $B$ & " + num2str(m_res{9}.B*12, '%0.3f') + " \\ \hline \end{tabular}"; 
writematrix(s5, "Output/LaTeX_Output/financial_frictions.txt", 'QuoteStrings', false);


s6 = "\begin{tabular}{lcc}  \hline \hline Moments & Model  & Data \\ " + ...
      " \hline Productivity growth & " + num2str(100*((1+m_res{11}.growth_rate)^12-1), '%0.2f') + "\% & 0.76\% \\ " + ...
      " Markups \\  " + ...
      " \hspace{0.2in} Mean & " + num2str(100*m_res{11}.markup, '%0.2f') + "\% & 19.41\% \\ " + ...
      " \hspace{0.2in} 50th percentile & " + num2str(100*(m_res{11}.p50_markup-1), '%0.2f') + "\% & 13.62\% \\ " + ...
      " \hspace{0.2in} 90th percentile & " + num2str(100*(m_res{11}.p90_markup-1), '%0.2f') + "\% & 42.65\% \\ " + ...
      " Profit volatility \\  \hspace{0.2in} All firms & " + num2str(100*m_res{11}.std_table.mean_fun1_profitgrowth_w(end), '%0.2f') + "\% & 35.72\% \\ " + ... 
      " \hspace{0.2in} Top profit quintile & " + num2str(100*m_res{11}.std_table.mean_fun1_profitgrowth_w(1), '%0.2f') + "\% & 19.48\% \\ " + ...
      " R\&D to sales \\  \hspace{0.2in} All firms & " + num2str(100*m_res{11}.std_table.mean_fun2_rdsales_w(end), '%0.2f') + "\% & 4.37\% \\ " + ...
      " \hspace{0.2in} Top profit quintile & " + num2str(100*m_res{11}.std_table.mean_fun2_rdsales_w(1), '%0.2f') + "\% & 2.64\% \\ " + ...
      " Innovation output \\  " + ...
      " \hspace{0.2in} Mean & " + num2str(100*m_res{11}.uncond_IO(1), '%0.2f') + "\% & 6.78\% \\ " + ...
      " \hspace{0.2in} 90th percentile & " + num2str(100*m_res{11}.uncond_IO(end-2), '%0.2f') + "\% & 19.54\% \\ " + ...
      " \hline \hline Parameters \\  \hline $B$ & " + num2str(12*m_res{11}.B, '%0.3f') + ...
      " \\  $\lambda$ & " + num2str(m_res{11}.lambda, '%0.3f') + ...
      " \\  $\eta$ & " + num2str(m_res{11}.eta(2)*12, '%0.3f') + ...
      " \\  \hline \end{tabular}"; 
writematrix(s6, "Output/LaTeX_Output/only_incremental.txt", 'QuoteStrings', false);

chi = 2.143390;
rk = 5.959426;
mtmp = benchmark_new();
mtmp = workhorse_nested_robust(mtmp);
rhotilde = rk - chi  - 100*((1 + mtmp.growth_rate)^12-1); 

chivec = repmat(sort([linspace(0, 0.003, 25), (1+chi/100)^(1/12)-1]), [1,2]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = 1:2
    for ii = 1:(length(chivec)/2) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = chivec(ii) + (1+rhotilde/100)^(1/12)-1;
        m_cell{count} = mtmp;
    end
end
gvec = zeros(1, length(m_cell)); 
rfvec = zeros(size(gvec)); 
markupvec = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec(ii) = 100*mtmp.markup; 
    rfvec(ii) = (1+rhotilde/100)^(1/12)-1 - chivec(ii) + mtmp.growth_rate; 
    rfvec(ii) = 100*((1+rfvec(ii))^12-1);
end



chivec = repmat(sort([linspace(8.46e-4, 0.004, 25), (1+chi/100)^(1/12)-1]), [1,2]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = 1:2
    for ii = 1:(length(chivec)/2) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = (1+chi/100)^(1/12)-1 +  chivec(ii); 
        m_cell{count} = mtmp;
    end
end
gvec_2 = zeros(1, length(m_cell)); 
rfvec_2 = zeros(size(gvec)); 
markupvec_2 = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec_2(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec_2(ii) = 100*mtmp.markup; 
    rfvec_2(ii) = -(1+chi/100)^(1/12)+1 + chivec(ii) + mtmp.growth_rate; 
    rfvec_2(ii) = 100*((1+rfvec_2(ii))^12-1);

end



close all
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=6;
y_width=3;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(1,2,1)
p1 = plot(rfvec(1:(length(chivec)/2)), gvec(1:(length(chivec)/2)), '-k', 'LineWidth', 2);
hold on 
p2 = plot(rfvec((length(chivec)/2+1):length(chivec)), gvec((length(chivec)/2+1):length(chivec)), '-r', ...
          'LineWidth', 2);
p3 = plot(rfvec_2(1:(length(chivec)/2)), gvec_2(1:(length(chivec)/2)), '--k', 'LineWidth', 2);
hold on 
p4 = plot(rfvec_2((length(chivec)/2+1):length(chivec)), gvec_2((length(chivec)/2+1):length(chivec)), '--r', ...
          'LineWidth', 2);
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('Growth, $g$', 'Interpreter', 'latex');
l1 = legend("No pos. R\&D, varying $\chi$", "No 3 qtr. profit, varying $\chi$", ...
            "No pos. R\&D, varying $\rho$", "No 3 qtr. profit, varying $\rho$", ...
            'Interpreter', 'latex', 'Location', 'north');
l1.FontSize = l1.FontSize*0.8;
set(l1, 'box', 'off');

subplot(1,2,2)
p1 = plot(rfvec(1:(length(chivec)/2)), markupvec(1:(length(chivec)/2)), '-k', 'LineWidth', 2);
hold on 
p2 = plot(rfvec((length(chivec)/2+1):length(chivec)), markupvec((length(chivec)/2+1):length(chivec)), '-r', ...
          'LineWidth', 2);
p3 = plot(rfvec_2(1:(length(chivec)/2)), markupvec_2(1:(length(chivec)/2)), '--k', 'LineWidth', 2);
hold on 
p4 = plot(rfvec_2((length(chivec)/2+1):length(chivec)), markupvec_2((length(chivec)/2+1):length(chivec)), '--r', ...
          'LineWidth', 2);
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('Mean net markup', 'Interpreter', 'latex');

saveas(gcf, "Output/Figures_Paper/moment_alt_fig.eps", 'epsc')
saveas(gcf, "Output/Figures_Paper/moment_alt_fig.png")


chi = 2.143390;
rk = 5.959426;
mtmp = benchmark_new();
mtmp = workhorse_nested_robust(mtmp);
rhotilde = rk - chi  - 100*((1+mtmp.growth_rate)^12-1); 
chivec = repmat(sort([linspace(0, 0.003, 25), (1+chi/100)^(1/12)-1]), [1,1]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = 3
    for ii = 1:(length(chivec)/1) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = chivec(ii) + (1+rhotilde/100)^(1/12)-1;
        m_cell{count} = mtmp;
    end
end
gvec = zeros(1, length(m_cell)); 
rfvec = zeros(size(gvec)); 
markupvec = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec(ii) = 100*(mtmp.markup_wt-1); 
    rfvec(ii) = (1+rhotilde/100)^(1/12)-1 - chivec(ii) + mtmp.growth_rate; 
    rfvec(ii) = 100*((1+rfvec(ii))^12-1);
end



chivec = repmat(sort([linspace(8.46e-4, 0.004, 25), (1+chi/100)^(1/12)-1]), [1,1]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = 3
    for ii = 1:(length(chivec)/1) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = (1+chi/100)^(1/12)-1 +  chivec(ii); 
        m_cell{count} = mtmp;
    end
end
gvec_2 = zeros(1, length(m_cell)); 
rfvec_2 = zeros(size(gvec)); 
markupvec_2 = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec_2(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec_2(ii) = 100*(mtmp.markup_wt-1); 
    rfvec_2(ii) = -(1+chi/100)^(1/12)+1 + chivec(ii) + mtmp.growth_rate; 
    rfvec_2(ii) = 100*((1+rfvec_2(ii))^12-1);

end

close all
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=3;
y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(2,1,1)
p1 = plot(rfvec(1:(length(chivec)/1)), gvec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), gvec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
hold on 
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('Growth, $g$', 'Interpreter', 'latex');
l1 = legend("Varying $\chi$", ...
            "Varying $\rho$", ...
            'Interpreter', 'latex', 'Location', 'north');
l1.FontSize = l1.FontSize*0.8;
set(l1, 'box', 'off');

subplot(2,1,2)
p1 = plot(rfvec(1:(length(chivec)/1)), markupvec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), markupvec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('Mean net markup', 'Interpreter', 'latex');
% ylim([17,22])

saveas(gcf, "Output/Figures_Paper/multistep_grho.eps", 'epsc')
saveas(gcf, "Output/Figures_Paper/multistep_grho.png")


chi = 2.143390;
rk = 5.959426;
mtmp = benchmark_new();
mtmp = workhorse_nested_robust(mtmp);
rhotilde = rk - chi  - 100*((1+mtmp.growth_rate)^12-1); 

chivec = repmat(sort([linspace(0, 0.003, 25), (1+chi/100)^(1/12)-1]), [1,1]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = 7
    for ii = 1:(length(chivec)/1) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = chivec(ii) + (1+rhotilde/100)^(1/12)-1;
        m_cell{count} = mtmp;
    end
end
gvec = zeros(1, length(m_cell)); 
rfvec = zeros(size(gvec)); 
markupvec = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec(ii) = 100*(mtmp.markup); 
    rfvec(ii) = (1+rhotilde/100)^(1/12)-1 - chivec(ii) + mtmp.growth_rate; 
    rfvec(ii) = 100*((1+rfvec(ii))^12-1);
end



chivec = repmat(sort([linspace(8.46e-4, 0.004, 25), (1+chi/100)^(1/12)-1]), [1,1]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = 7
    for ii = 1:(length(chivec)/1) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = (1+chi/100)^(1/12)-1 +  chivec(ii); 
        m_cell{count} = mtmp;
    end
end
gvec_2 = zeros(1, length(m_cell)); 
rfvec_2 = zeros(size(gvec)); 
markupvec_2 = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec_2(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec_2(ii) = 100*(mtmp.markup); 
    rfvec_2(ii) = -(1+chi/100)^(1/12)+1 + chivec(ii) + mtmp.growth_rate; 
    rfvec_2(ii) = 100*((1+rfvec_2(ii))^12-1);

end

close all
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=3;
y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(2,1,1)
p1 = plot(rfvec(1:(length(chivec)/1)), gvec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), gvec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
hold on 
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('Growth, $g$', 'Interpreter', 'latex');
l1 = legend("Varying $\chi$", ...
            "Varying $\rho$", ...
            'Interpreter', 'latex', 'Location', 'north');
l1.FontSize = l1.FontSize*0.8;
set(l1, 'box', 'off');

subplot(2,1,2)
p1 = plot(rfvec(1:(length(chivec)/1)), markupvec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), markupvec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('Mean net markup', 'Interpreter', 'latex');
% ylim([17,22])s

saveas(gcf, "Output/Figures_Paper/highmarkup_grho.eps", 'epsc')
saveas(gcf, "Output/Figures_Paper/highmarkup_grho.png")





chi = 2.143390;
rk = 5.959426;
mtmp = benchmark_new();
mtmp = workhorse_nested_robust(mtmp);
rhotilde = rk - chi  - 100*((1+mtmp.growth_rate)^12-1); 

chivec = repmat(sort([linspace(0, 0.003, 25), (1+chi/100)^(1/12)-1]), [1,4]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = [4, 8, 5, 6]
    for ii = 1:(length(chivec)/4) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = chivec(ii) + (1+rhotilde/100)^(1/12)-1;
        m_cell{count} = mtmp;
    end
end
gvec = zeros(1, length(m_cell)); 
rfvec = zeros(size(gvec)); 
markupvec = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec(ii) = 100*(mtmp.markup); 
    rfvec(ii) = (1+rhotilde/100)^(1/12)-1 - chivec(ii) + mtmp.growth_rate; 
    rfvec(ii) = 100*((1+rfvec(ii))^12-1);
end



chivec = repmat(sort([linspace(8.46e-4, 0.004, 25), (1+chi/100)^(1/12)-1]), [1,4]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = [4, 8, 5, 6]
    for ii = 1:(length(chivec)/4) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = (1+chi/100)^(1/12)-1 +  chivec(ii); 
        m_cell{count} = mtmp;
    end
end
gvec_2 = zeros(1, length(m_cell)); 
rfvec_2 = zeros(size(gvec)); 
markupvec_2 = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec_2(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec_2(ii) = 100*(mtmp.markup); 
    rfvec_2(ii) = -(1+chi/100)^(1/12)+1 + chivec(ii) + mtmp.growth_rate; 
    rfvec_2(ii) = 100*((1+rfvec_2(ii))^12-1);

end


close all 
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=6.5;
y_width=4.5;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(2,2,1)
p1 = plot(rfvec(1:(length(chivec)/4)), gvec(1:(length(chivec)/4)), '-k', 'LineWidth', 2);
hold on 
p2 = plot(rfvec((length(chivec)/4 + 1):(2*length(chivec)/4)), gvec((length(chivec)/4+1):(2*length(chivec)/4)), '--r', 'LineWidth', 2);
p3 = plot(rfvec((2*length(chivec)/4 + 1):(3*length(chivec)/4)), gvec((2*length(chivec)/4+1):(3*length(chivec)/4)), ...
          'LineStyle', ':', 'color', [0, 0, 1], 'LineWidth', 2);
p4 = plot(rfvec((3*length(chivec)/4+1):length(chivec)), gvec((3*length(chivec)/4 +1):length(chivec)), 'LineStyle', '-.', ...
          'color', [0, 0.5, 0], 'LineWidth', 2);
title({'Varying $\chi$'; ''}, 'Interpreter', 'latex')
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex')
ylabel('Growth rate, $g$', 'Interpreter', 'latex')

subplot(2,2,2)
p1 = plot(rfvec_2(1:(length(chivec)/4)), gvec_2(1:(length(chivec)/4)), '-k', 'LineWidth', 2);
hold on 
p2 = plot(rfvec_2((length(chivec)/4 + 1):(2*length(chivec)/4)), gvec_2((length(chivec)/4+1):(2*length(chivec)/4)), '--r', 'LineWidth', 2);
p3 = plot(rfvec_2((2*length(chivec)/4 + 1):(3*length(chivec)/4)), gvec_2((2*length(chivec)/4+1):(3*length(chivec)/4)), ...
          'LineStyle', ':', 'color', [0, 0, 1], 'LineWidth', 2);
p4 = plot(rfvec_2((3*length(chivec)/4+1):length(chivec)), gvec_2((3*length(chivec)/4 +1):length(chivec)), 'LineStyle', '-.', ...
          'color', [0, 0.5, 0], 'LineWidth', 2);
title({'Varying $\rho$'; ''}, 'Interpreter', 'latex')
ylabel('Growth rate, $g$', 'Interpreter', 'latex')
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex')
l1 = legend('$\gamma = 0.35$', '$\eta = 0.03$', '$\varphi = 0.50$', 'Inelastic labor', ...
            'Interpreter', 'latex', 'Location', 'northeast');
set(l1, 'box', 'off');


subplot(2,2,3)
p1 = plot(rfvec(1:(length(chivec)/4)), markupvec(1:(length(chivec)/4)), '-k', 'LineWidth', 2);
hold on 
p2 = plot(rfvec((length(chivec)/4 + 1):(2*length(chivec)/4)), markupvec((length(chivec)/4+1):(2*length(chivec)/4)), '--r', 'LineWidth', 2);
p3 = plot(rfvec((2*length(chivec)/4 + 1):(3*length(chivec)/4)), markupvec((2*length(chivec)/4+1):(3*length(chivec)/4)), ...
          'LineStyle', ':', 'color', [0, 0, 1], 'LineWidth', 2);
p4 = plot(rfvec((3*length(chivec)/4+1):length(chivec)), markupvec((3*length(chivec)/4 +1):length(chivec)), 'LineStyle', '-.', ...
          'color', [0, 0.5, 0], 'LineWidth', 2);
ylabel('Mean net markup', 'Interpreter', 'latex')
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex')


subplot(2,2,4)
p1 = plot(rfvec_2(1:(length(chivec)/4)), markupvec_2(1:(length(chivec)/4)), '-k', 'LineWidth', 2);
hold on 
p2 = plot(rfvec_2((length(chivec)/4 + 1):(2*length(chivec)/4)), markupvec_2((length(chivec)/4+1):(2*length(chivec)/4)), '--r', 'LineWidth', 2);
p3 = plot(rfvec_2((2*length(chivec)/4 + 1):(3*length(chivec)/4)), markupvec_2((2*length(chivec)/4+1):(3*length(chivec)/4)), ...
          'LineStyle', ':', 'color', [0, 0, 1], 'LineWidth', 2);
p4 = plot(rfvec_2((3*length(chivec)/4+1):length(chivec)), markupvec_2((3*length(chivec)/4 +1):length(chivec)), 'LineStyle', '-.', ...
          'color', [0, 0.5, 0], 'LineWidth', 2);
ylabel('Mean net markup', 'Interpreter', 'latex')
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex')


saveas(gcf, "Output/Figures_Paper/robust_2x2.eps", 'epsc')
saveas(gcf, "Output/Figures_Paper/robust_2x2.png")

chi = 2.143390;
rk = 5.959426;
mtmp = benchmark_new();
mtmp = workhorse_nested_robust(mtmp);
rhotilde = rk - chi  - 100*((1+mtmp.growth_rate)^12-1); 

chivec = repmat(sort([linspace(0, 0.003, 25), (1+chi/100)^(1/12)-1]), [1,1]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = 10
    for ii = 1:(length(chivec)/1) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = chivec(ii) + (1+rhotilde/100)^(1/12)-1;
        m_cell{count} = mtmp;
    end
end
gvec = zeros(1, length(m_cell)); 
rfvec = zeros(size(gvec)); 
markupvec = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec(ii) = 100*(mtmp.markup_wt-1); 
    rfvec(ii) = (1+rhotilde/100)^(1/12)-1 - chivec(ii) + mtmp.growth_rate; 
    rfvec(ii) = 100*((1+rfvec(ii))^12-1);
end



chivec = repmat(sort([linspace(8.46e-4, 0.004, 25), (1+chi/100)^(1/12)-1]), [1,1]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = 10
    for ii = 1:(length(chivec)/1) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = (1+chi/100)^(1/12)-1 +  chivec(ii); 
        m_cell{count} = mtmp;
    end
end
gvec_2 = zeros(1, length(m_cell)); 
rfvec_2 = zeros(size(gvec)); 
markupvec_2 = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec_2(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec_2(ii) = 100*(mtmp.markup_wt-1); 
    rfvec_2(ii) = -(1+chi/100)^(1/12)+1 + chivec(ii) + mtmp.growth_rate; 
    rfvec_2(ii) = 100*((1+rfvec_2(ii))^12-1);

end


close all
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=3;
y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(2,1,1)
p1 = plot(rfvec(1:(length(chivec)/1)), gvec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), gvec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
hold on 
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('Growth, $g$', 'Interpreter', 'latex');
l1 = legend("Varying $\chi$", ...
            "Varying $\rho$", ...
            'Interpreter', 'latex', 'Location', 'north');
l1.FontSize = l1.FontSize*0.8;
set(l1, 'box', 'off');

subplot(2,1,2)
p1 = plot(rfvec(1:(length(chivec)/1)), markupvec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), markupvec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('Mean net markup', 'Interpreter', 'latex');


saveas(gcf, "Output/Figures_Paper/kappa12_grho.eps", 'epsc')
saveas(gcf, "Output/Figures_Paper/kappa12_grho.png")


chi = 2.143390;
rk = 5.959426;
mtmp = benchmark_new();
mtmp = workhorse_nested_robust(mtmp);
rhotilde = rk - chi  - 100*((1+mtmp.growth_rate)^12-1); 

chivec = repmat(sort([linspace(0, 0.003, 25), (1+chi/100)^(1/12)-1]), [1,1]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = 11
    for ii = 1:(length(chivec)/1) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = chivec(ii) + (1+rhotilde/100)^(1/12)-1;
        m_cell{count} = mtmp;
    end
end
gvec = zeros(1, length(m_cell)); 
rfvec = zeros(size(gvec)); 
markupvec = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec(ii) = 100*(mtmp.markup); 
    rfvec(ii) = (1+rhotilde/100)^(1/12)-1 - chivec(ii) + mtmp.growth_rate; 
    rfvec(ii) = 100*((1+rfvec(ii))^12-1);
end



chivec = repmat(sort([linspace(8.46e-4, 0.004, 25), (1+chi/100)^(1/12)-1]), [1,1]);
m_cell = cell(1, length(chivec));
count = 0;
for jj = 11
    for ii = 1:(length(chivec)/1) 
        count = count + 1; 
        mtmp = mvec{jj}();
        mtmp.rho = (1+chi/100)^(1/12)-1 +  chivec(ii); 
        m_cell{count} = mtmp;
    end
end
gvec_2 = zeros(1, length(m_cell)); 
rfvec_2 = zeros(size(gvec)); 
markupvec_2 = zeros(size(gvec)); 
parfor (ii = 1:length(m_cell))

    mtmp = workhorse_nested_robust(m_cell{ii}); 
    gvec_2(ii) = 100*((1+mtmp.growth_rate)^12-1);
    markupvec_2(ii) = 100*(mtmp.markup); 
    rfvec_2(ii) = -(1+chi/100)^(1/12)+1 + chivec(ii) + mtmp.growth_rate; 
    rfvec_2(ii) = 100*((1+rfvec_2(ii))^12-1);

end


close all
figure; 

set(gcf, 'PaperUnits', 'inches');
x_width=3;
y_width=4;
set(gcf, 'PaperPosition', [0 0 x_width y_width]); %


subplot(2,1,1)
p1 = plot(rfvec(1:(length(chivec)/1)), gvec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), gvec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
hold on 
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('Growth, $g$', 'Interpreter', 'latex');
ylim([0.6, 1.2])
l1 = legend("Varying $\chi$", ...
            "Varying $\rho$", ...
            'Interpreter', 'latex', 'Location', 'north');
l1.FontSize = l1.FontSize*0.8;
set(l1, 'box', 'off');

subplot(2,1,2)
p1 = plot(rfvec(1:(length(chivec)/1)), markupvec(1:(length(chivec)/1)), '-k', 'LineWidth', 2);
hold on 
p3 = plot(rfvec_2(1:(length(chivec)/1)), markupvec_2(1:(length(chivec)/1)), '--k', 'LineWidth', 2);
xlabel('Risk-free rate, $r^f$', 'Interpreter', 'latex');
ylabel('Mean net markup', 'Interpreter', 'latex');
% ylim([17,22])s

saveas(gcf, "Output/Figures_Paper/onlyincremental_grho.eps", 'epsc')
saveas(gcf, "Output/Figures_Paper/onlyincremental_grho.png")


end