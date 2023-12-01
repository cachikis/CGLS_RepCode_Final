%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: moment_sensitivity_new.m
% Author: Craig A. Chikis
% Date: 08/07/2023
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [] = moment_sensitivity_new()

	   
	m = benchmark_new();

	B_vec = m.B*[0.05,.1,.25,.5,.85,.9,.95,1,1.05,1.1, 1.15,1.5,1.75,1.9,2.5];
	lambda_vec = 1 + (m.lambda-1)*[0.05,.1,.25,.5,.85,.9,.95,1,1.05,1.1, 1.15,1.5,1.75,1.9,2.5];
	phi_vec = min(m.phi_wt*[0.05,.1,.25,.5,.85,.9,.95,1,1.05,1.1, 1.15,1.5,1.75,1.9,2.5],1);
	gamma_vec = sort([linspace(0.05, 0.75, length(phi_vec)-1), m.gamma]);

	grid_points = (length(B_vec)-1)/2;


	m_collect = cell(4, length(B_vec)); 
	for (oi = 1:size(m_collect, 1)) 
		for (ii = 1:size(B_vec, 2))

			m = benchmark_new(); 
			if (oi == 1)
				m.B = B_vec(ii);
			elseif (oi == 2)
				m.lambda = lambda_vec(ii);
			elseif (oi == 3)
				m.phi_wt = phi_vec(ii);
				[m.prob, ~, ~] = F_mat_main(m.phi_wt, m.phi_wt_tilde, m.phi_wt_exog, m.l, m.l_tilde, m.s_max);
			elseif (oi == 4)
				m.gamma = gamma_vec(ii);
			end

			m = workhorse_nested_robust(m);
			try 
				m = sim_wrap(m);
				m = CSTAT_sim(m); 
				m = innovation_output(m); 
                m = FHK_sim(m); 
			catch
			end
			[m.p50_approx, ~, ~] = markup_quant(m.kappa, m.lambda, m.s_max, m.firm_distributions, m.nu_s, 0.5);
			[m.p90_approx, m.markup_array, ~] = markup_quant(m.kappa, m.lambda, m.s_max, m.firm_distributions, m.nu_s, 0.9);
	
			m_collect{oi, ii} = m; 
		end


	end

	growth_vec = nan(size(m_collect,1), size(m_collect,2));

	markup_vec = nan(size(growth_vec));
	markup_p50_vec = nan(size(growth_vec));
	markup_p75_vec = nan(size(growth_vec));
	markup_p90_vec = nan(size(growth_vec));

	within_vec = nan(size(growth_vec)); 

	innov_prof_vol_vec = nan(size(growth_vec));
	innov_prof_vol_quint_vec = nan(size(growth_vec));

	RD_sales_median_vec = nan(size(growth_vec));
	RD_sales_quint_vec = nan(size(growth_vec));

	innov_mean_vec = nan(size(growth_vec));
	innov_p90_vec = nan(size(growth_vec));
	for (oi = 1:size(m_collect,1))
		for (ii = 1:size(m_collect,2))
			growth_vec(oi, ii) = 100*(power(1+m_collect{oi, ii}.growth_rate, 12)-1); 

			markup_vec(oi,ii) = 100*m_collect{oi, ii}.markup; 
			markup_p90_vec(oi,ii) = 100*m_collect{oi, ii}.p90_approx;

			try 
				innov_prof_vol_vec(oi,ii) = 100*m_collect{oi,ii}.std_table.mean_fun1_profitgrowth_w(end);
				innov_prof_vol_quint_vec(oi,ii) = 100*m_collect{oi,ii}.std_table.mean_fun1_profitgrowth_w(1);

				RD_sales_median_vec(oi,ii) = 100*m_collect{oi,ii}.std_table.mean_fun2_rdsales_w(end);
				RD_sales_quint_vec(oi,ii) = 100*m_collect{oi,ii}.std_table.mean_fun2_rdsales_w(1);

				within_vec(oi,ii) = 100*m_collect{oi,ii}.WITHIN_out;

				innov_mean_vec(oi, ii) = 100*m_collect{oi,ii}.uncond_IO(1);
				innov_p90_vec(oi, ii) = 100*m_collect{oi,ii}.uncond_IO(9);
			catch
			end
		end
	end

	m = benchmark_new();
	param_store = {growth_vec, markup_vec, ...
				   innov_prof_vol_vec, innov_prof_vol_quint_vec, RD_sales_median_vec, RD_sales_quint_vec, ...
				   innov_mean_vec, innov_p90_vec, within_vec};

	idx_choose = round(linspace(1, size(param_store{1},2), 7));
	idx_choose2 = [2, 8, 14]; 
	

	if (~ismember(m.phi_wt, phi_vec(idx_choose2)))
		idx_choose2 = sort([idx_choose2, 8]);
	end
	if (~ismember(m.phi_wt, phi_vec(idx_choose)))
		idx_choose = sort([idx_choose, 8]);
	end


	change_cell = cell(1, length(param_store));
	change_cell2 = cell(size(change_cell)); 
	mod_store = cell(size(param_store{1},1), length(param_store));
	mod_store2 = cell(size(mod_store)); 
	
	for (oi = 1:length(param_store))
		mom_change = zeros(size(param_store{oi},1), 1);
		mom_change2 = zeros(size(mom_change)); 
		for (ii = 1:size(param_store{oi},1))
			if (ii == 1)
				x = B_vec;
				param_fix = m.B;
				param_inc = m.B*1.05;
			elseif (ii == 2)
				x = lambda_vec;
				param_fix = m.lambda;
				param_inc = 1 + (m.lambda-1)*1.05;
			elseif (ii == 3)
				x = phi_vec;
				param_fix = m.phi_wt;
				param_inc = m.phi_wt*1.05;
			elseif (ii == 4)
				x = gamma_vec;
				param_fix = m.gamma;
				param_inc = m.gamma*1.05;
			end

			y = param_store{oi}(ii, :);
			y1 = y(idx_choose); 
			x1 = x(idx_choose); 
			y2 = y(idx_choose2);
			x2 = x(idx_choose2); 
			remove1 = find(isnan(y1) | isnan(x1)); 
			remove2 = find(isnan(y2) | isnan(x2)); 

			x1(remove1) = [];
			y1(remove1) = [];
			x2(remove2) = [];
			y2(remove2) = []; 
		
			model = spline(x1, y1);
			model2 = spline(x2, y2);  
			mod_store{ii, oi} = model;
			mod_store2{ii, oi} = model2; 

			mom_change(ii) = round(100*(ppval(model, param_inc)/ppval(model, param_fix) - 1), 3, 'significant');
			mom_change2(ii) = round(100*(ppval(model2, param_inc)/ppval(model2, param_fix) - 1), 3, 'significant');
		end

		change_cell{oi} = mom_change;
		change_cell2{oi} = mom_change2; 
	end

	growth_change = change_cell{1}';
	markup_change = change_cell{2}'; 
	innov_prof_vol_change = change_cell{3}';
	innov_prof_vol_quint_change = change_cell{4}';
	RD_sales_median_change = change_cell{5}';
	RD_sales_quint_change = change_cell{6}';
	innov_mean_change = change_cell{7}';
	innov_p90_change_s = change_cell2{8}';
	


	s1 = strcat("\begin{tabular}{lcccc}  \hline \hline &  \multicolumn{3}{c}{Change in ...} & \\ \cmidrule(lr){2-4} ", ...
				" \multicolumn{1}{l}{Moments} & \multicolumn{1}{c}{$\phi$} & \multicolumn{1}{c}{$\lambda$} & ", ...
				" \multicolumn{1}{c}{$B$} & \multicolumn{1}{l}{Initial level} \\ \hline ", ...
				"\shortstack{Productivity growth} & ", num2str(growth_change(1,3), '%0.2f'), "\%", " & ", num2str(growth_change(1,2), '%0.2f'), "\%", " & ", ...
				num2str(growth_change(1,1), '%0.2f'),  "\%", " & ", num2str(growth_vec(1,8), '%0.2f'), "\%", " \\ ", ...
							...
				"  Mean net markup & ", num2str(markup_change(1,3), '%0.2f'), "\%", " & ", num2str(markup_change(1,2), '%0.2f'),"\%", " & ", ...
				num2str(markup_change(1,1), '%0.2f'), "\%", " & ", num2str(markup_vec(1,8),'%0.2f'), "\%", " \\ ", ...
							...
				"Profit volatility & & & & \\ ", ...
				"\hspace{0.2in} All firms & ", num2str(innov_prof_vol_change(1,3), '%0.2f'), "\%"," & ", num2str(innov_prof_vol_change(1,2),'%0.2f'), "\%", " & ", ...
				num2str(innov_prof_vol_change(1,1), '%0.2f'),  "\%", " & ", num2str(innov_prof_vol_vec(1,8), '%0.2f'), "\%", " \\ ", ...
				"\hspace{0.2in} Top profit quintile & ", num2str(innov_prof_vol_quint_change(1,3), '%0.2f'), "\%"," & ", num2str(innov_prof_vol_quint_change(1,2), '%0.2f'),"\%", " & ", ...
				num2str(innov_prof_vol_quint_change(1,1), '%0.2f'),  "\%", " & ", num2str(innov_prof_vol_quint_vec(1,8), '%0.2f'), "\%", " \\ ", ...
				"R\&D to sales & & & & \\ ", ...
				"\hspace{0.2in} All firms & ", num2str(RD_sales_median_change(1,3), '%0.2f'), "\%"," & ", num2str(RD_sales_median_change(1,2), '%0.2f'), "\%"," & ", ...
				num2str(RD_sales_median_change(1,1), '%0.2f'),  "\%", " & ", num2str(RD_sales_median_vec(1,8), '%0.2f'), "\%", " \\ ", ...
				"\hspace{0.2in} Top profit quintile & ", num2str(RD_sales_quint_change(1,3), '%0.2f'),"\%", " & ", num2str(RD_sales_quint_change(1,2), '%0.2f'),"\%", " & ", ...
				num2str(RD_sales_quint_change(1,1), '%0.2f'),  "\%", " & ", num2str(RD_sales_quint_vec(1,8),'%0.2f'), "\%", " \\ ", ...
				"Innovation output & & & & \\ ", ...
				"\hspace{0.2in} Mean & ", num2str(innov_mean_change(1,3), '%0.2f'), "\%"," & ", num2str(innov_mean_change(1,2),'%0.2f'), "\%", " & ", ...
				num2str(innov_mean_change(1,1), '%0.2f'),  "\%", " & ", num2str(innov_mean_vec(1,8), '%0.2f'), "\%", " \\ ", ...
				"\hspace{0.2in} 90th percentile & ", num2str(innov_p90_change_s(1,3), '%0.2f'), "\%"," & ", num2str(innov_p90_change_s(1,2), '%0.2f'),"\%", " & ", ...
				num2str(innov_p90_change_s(1,1), '%0.2f'),  "\%", " & ", num2str(innov_p90_vec(1,8),'%0.2f'), "\%", " ", ...
				" \\ \hline ", " \end{tabular}");

	writematrix(s1, "Output/LaTeX_Output/sensitivity_table.txt", ...
			    'QuoteStrings', false);


	rng(2022+08+30, 'Threefry')
    draw_Lerner = betarnd(1.36, 8, 1e6, 1);
    draw_markup = 1./(1-draw_Lerner) - 1; 

    target_filter_vec = ["three_dropm1"]; 
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

	rdsales_targets = readtable(strcat("Output/Store_Data/rdsales_t_targets_", target_filter_vec(1), ".csv")); 
    profvol_targets = readtable(strcat("Output/Store_Data/profitvol_t_", target_filter_vec(1), ".csv")); 
	uncond_IO_targets = readmatrix(strcat("Output/Store_Data/kogan_targets.csv"));
    markup_targets = [mean(draw_markup), prctile(draw_markup, [50, 90])]; 
    growth_targ = 1.033; 

	clear draw_markup; 
	clear draw_Lerner;
					
   
	multiples = [.1,.25,.5,.85,.95,1,1.05, 1.15,1.5,1.75,1.9,2.5];
	plot_cell = {[5:9], [1:length(B_vec)]};
	print_fig = ["local", "global"];

	ii = 2;
	plot_vec = plot_cell{ii}(2:(end-1));

	close all
	figure; 

	B_fix = m.B;

	m = benchmark_new();

	targets = [growth_targ, markup_targets(1)*100, rdsales_targets.p50(end)*100, profvol_targets.sd(end)*100, ...
				uncond_IO_targets(1)]; 

	set(gcf, 'PaperUnits', 'inches');
	x_width=6;
	y_width=7.5;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %

	subplot(5,3,1)
	plot(phi_vec(plot_vec), growth_vec(3,plot_vec), '-b', 'LineWidth', 2);
	hold on 
	plot(linspace(m.phi_wt - 0.25*(max(phi_vec) - min(phi_vec)), ...
				  m.phi_wt + 0.25*(max(phi_vec) - min(phi_vec)), 10), ...
		 ones(1,10)*targets(1), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	ylabel({'Growth, $g$'; ''}, 'Interpreter', 'latex');
	xline(m.phi_wt, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([m.phi_wt - .25, m.phi_wt + .25]);
	ylim([0,2.5])
	xticks([0.1,0.2,0.3,0.4,0.5])
	yticks([0,0.5,1,1.5,2,2.5])


	subplot(5,3,2)
	plot(lambda_vec(plot_vec), growth_vec(2,plot_vec), '--r', 'LineWidth', 2);
	hold on 
	plot(linspace(m.lambda - 0.25*(m.lambda+0.015 - m.lambda + 0.015), ...
				  m.lambda + 0.25*(m.lambda+0.015 - m.lambda + 0.015), 10), ...
		 ones(1,10)*targets(1), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(m.lambda, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([m.lambda - .0155, m.lambda + .0155]);
	ylim([0,2.5])
	xticks([1.01,1.02,1.03,1.04])
	yticks([0,0.5,1,1.5,2,2.5])



	subplot(5,3,3)
	plot(12*(B_vec(plot_vec)*(1-0)^(-m.gamma)), growth_vec(1,plot_vec), ':k', 'LineWidth', 2);
	hold on 
	plot(linspace(12*(m.B*(1-0)^(-m.gamma)) - 0.25*((m.B*(1-0)^(-m.gamma))*12 + 2 - ((m.B*(1-0)^(-m.gamma))*12 - 2)), ...
				  12*(m.B*(1-0)^(-m.gamma)) + 0.25*((m.B*(1-0)^(-m.gamma))*12 + 2 - ((m.B*(1-0)^(-m.gamma))*12 - 2)), 10), ...
		 ones(1,10)*targets(1), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(B_fix*12, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([0, 2.5]);
	ylim([0,2.5])
	yticks([0,0.5,1,1.5,2,2.5])
	xticks([0,0.5,1,1.5,2,2.5])


	subplot(5,3,4)
	plot(phi_vec(plot_vec), markup_vec(3,plot_vec), '-b', 'LineWidth', 2);
	hold on 
	plot(linspace(m.phi_wt - 0.25*(max(phi_vec) - min(phi_vec)), ...
				  m.phi_wt + 0.25*(max(phi_vec) - min(phi_vec)), 10), ...
		 ones(1,10)*targets(2), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	ylabel({'Mean net markup'; ''}, 'Interpreter', 'latex');
	ylim([0, 400]);
	yticks([0, 200, 400]);
	xline(m.phi_wt, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([m.phi_wt - .25, m.phi_wt + .25]);
	xticks([0.1,0.2,0.3,0.4,0.5])

	subplot(5,3,5)
	plot(lambda_vec(plot_vec), markup_vec(2,plot_vec), '--r', 'LineWidth', 2);
	hold on 
	plot(linspace(m.lambda - 0.25*(m.lambda+0.015 - m.lambda + 0.015), ...
				  m.lambda + 0.25*(m.lambda+0.015 - m.lambda + 0.015), 10), ...
		 ones(1,10)*targets(2), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	ylim([0, 100]);
	yticks([0,50,100]);
	xline(m.lambda, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([m.lambda - .0155, m.lambda + .0155]);
	xticks([1.01,1.02,1.03,1.04])


	subplot(5,3,6)
	plot(12*(B_vec(plot_vec)*(1-0)^(-m.gamma)), markup_vec(1,plot_vec), ':k', 'LineWidth', 2);
	hold on 
	plot(linspace(12*(m.B*(1-0)^(-m.gamma)) - 0.25*((m.B*(1-0)^(-m.gamma))*12 + 2 - ((m.B*(1-0)^(-m.gamma))*12 - 2)), ...
				  12*(m.B*(1-0)^(-m.gamma)) + 0.25*((m.B*(1-0)^(-m.gamma))*12 + 2 - ((m.B*(1-0)^(-m.gamma))*12 - 2)), 10), ...
		 ones(1,10)*targets(2), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(B_fix*12, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([0, 2.5]);
	ylim([0, 100]);
	yticks([0,50,100]);
	xticks([0,0.5,1,1.5,2,2.5])

	subplot(5,3,7)
	plot(phi_vec(plot_vec), RD_sales_median_vec(3,plot_vec), '-b', 'LineWidth', 2);
	hold on 
	plot(linspace(m.phi_wt - 0.25*(max(phi_vec) - min(phi_vec)), ...
				  m.phi_wt + 0.25*(max(phi_vec) - min(phi_vec)), 10), ...
		 ones(1,10)*targets(3), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	ylabel({'R\&D to sales'; '(median, all firms)'}, 'Interpreter', 'latex');
	xline(m.phi_wt, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([m.phi_wt - .25, m.phi_wt + .25]);
	ylim([0,8])
	xticks([0.1,0.2,0.3,0.4,0.5])
	yticks([0,4,8])

	subplot(5,3,8)
	plot(lambda_vec(plot_vec),RD_sales_median_vec(2,plot_vec),'--r', 'LineWidth', 2);
	hold on 
	plot(linspace(m.lambda - 0.25*(m.lambda+0.015 - m.lambda + 0.015), ...
				  m.lambda + 0.25*(m.lambda+0.015 - m.lambda + 0.015), 10), ...
		 ones(1,10)*targets(3), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(m.lambda, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([m.lambda - .0155, m.lambda + .0155]);
	xticks([1.01,1.02,1.03,1.04])
	yticks([0,4,8])
	ylim([0,8])


	subplot(5,3,9)
	plot((B_vec(plot_vec)*(1-0)^(-m.gamma))*12, RD_sales_median_vec(1,plot_vec), ':k', 'LineWidth', 2);
	hold on 
	plot(linspace(12*(m.B*(1-0)^(-m.gamma)) - 0.25*((m.B*(1-0)^(-m.gamma))*12 + 2 - ((m.B*(1-0)^(-m.gamma))*12 - 2)), ...
				  12*(m.B*(1-0)^(-m.gamma)) + 0.25*((m.B*(1-0)^(-m.gamma))*12 + 2 - ((m.B*(1-0)^(-m.gamma))*12 - 2)), 10), ...
		 ones(1,10)*targets(3),...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(B_fix*12, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([0, 2.5]);
	ylim([0,8])
	yticks([0,4,8])
	xticks([0,0.5,1,1.5,2,2.5])


	subplot(5,3,10)
	plot(phi_vec(plot_vec), innov_prof_vol_vec(3,plot_vec), '-b', 'LineWidth', 2);
	hold on 
	plot(linspace(m.phi_wt - 0.25*(max(phi_vec) - min(phi_vec)), ...
				  m.phi_wt + 0.25*(max(phi_vec) - min(phi_vec)), 10), ...
		 ones(1,10)*targets(4), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	ylabel({'Profit volatility'; '(all firms)'}, 'Interpreter', 'latex');
	xline(m.phi_wt, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([m.phi_wt - .25, m.phi_wt + .25]);
	ylim([0,50])
	yticks([0,25,50])
	xticks([0.1,0.2,0.3,0.4,0.5])

	subplot(5,3,11)
	plot(lambda_vec(plot_vec),innov_prof_vol_vec(2,plot_vec),'--r', 'LineWidth', 2);
	hold on 
	plot(linspace(m.lambda - 0.25*(m.lambda+0.015 - m.lambda + 0.015), ...
				  m.lambda + 0.25*(m.lambda+0.015 - m.lambda + 0.015), 10), ...
		 ones(1,10)*targets(4), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(m.lambda, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([m.lambda - .0155, m.lambda + .0155]);
	ylim([0,50])
	yticks([0,25,50])
	xticks([1.01,1.02,1.03,1.04])

	subplot(5,3,12)
	plot((B_vec(plot_vec)*(1-0)^(-m.gamma))*12, innov_prof_vol_vec(1,plot_vec), ':k', 'LineWidth', 2);
	hold on 
	plot(linspace(12*(m.B*(1-0)^(-m.gamma)) - 0.25*((m.B*(1-0)^(-m.gamma))*12 + 2 - ((m.B*(1-0)^(-m.gamma))*12 - 2)), ...
				  12*(m.B*(1-0)^(-m.gamma)) + 0.25*((m.B*(1-0)^(-m.gamma))*12 + 2 - ((m.B*(1-0)^(-m.gamma))*12 - 2)), 10), ...
		 ones(1,10)*targets(4), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(B_fix*12, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([0, 2.5]);
	ylim([0,50])
	yticks([0,25,50])
	xticks([0,0.5,1,1.5,2,2.5])

	phi_fix = m.phi_wt;
	phi_wt = phi_fix; 
	lambda_fix = m.lambda; 
	lambda = lambda_fix; 
	B_fix = m.B;
	B = B_fix; 
	gamma_fix = m.gamma;
	gamma = gamma_fix;  
	subsidy = m.subsidy; 

	subplot(5,3,13)
	plot(phi_vec(plot_vec), innov_mean_vec(3,plot_vec), '-b', 'LineWidth', 2);
	hold on 
	plot(linspace(m.phi_wt - 0.25*(max(phi_vec) - min(phi_vec)), ...
				  m.phi_wt + 0.25*(max(phi_vec) - min(phi_vec)), 10), ...
		 ones(1,10)*targets(5), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	ylabel({'Innovation output'; '(mean)'}, 'Interpreter', 'latex');
	xline(phi_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([phi_fix - .25, phi_fix + .25]);
	ylim([0,12])
	yticks([0,4,8,12])
	xticks([0.1,0.2,0.3,0.4,0.5])
	xlabel('$\phi$', 'Interpreter', 'latex')

	subplot(5,3,14)
	plot(lambda_vec(plot_vec),innov_mean_vec(2,plot_vec),'--r', 'LineWidth', 2);
	hold on 
	plot(linspace(lambda - 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), ...
				  lambda + 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), 10), ...
		 ones(1,10)*targets(5), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(lambda_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([lambda_fix - .0155, lambda_fix + .0155]);
	ylim([0,12])
	yticks([0,4,8,12])
	xticks([1.01,1.02,1.03,1.04])
	xlabel('$\lambda$', 'Interpreter', 'latex')

	subplot(5,3,15)
	plot((B_vec(plot_vec)*(1-0)^(-gamma))*12, innov_mean_vec(1,plot_vec), ':k', 'LineWidth', 2);
	hold on 
	plot(linspace(12*(B*(1-0)^(-gamma)) - 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), ...
				  12*(B*(1-0)^(-gamma)) + 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), 10), ...
		 ones(1,10)*targets(5), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(B_fix*12, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([0, 2.5]);
	ylim([0,12])
	yticks([0,4,8,12])
	xticks([0,.5,1,1.5,2,2.5])
	xlabel('$B$', 'Interpreter', 'latex')

	



	saveas(gcf, strcat("Output/Figures_Paper/fourmom_text.eps"), 'epsc');
	saveas(gcf, strcat("Output/Figures_Paper/fourmom_text.png"));




	
	close all
	figure; 

	benchmark_new();

	targets = [nan, markup_targets(3)*100, rdsales_targets.p50(1)*100, profvol_targets.sd(1)*100, FHK_targ(1), ...
				uncond_IO_targets(9)];
		
	set(gcf, 'PaperUnits', 'inches');
	x_width=6;
	y_width=7.5;
	set(gcf, 'PaperPosition', [0 0 x_width y_width]); %



	subplot(5,4,1)
	plot(phi_vec(plot_vec), markup_p90_vec(3,plot_vec)-100, '-b', 'LineWidth', 2);
	hold on 
	plot(linspace(phi_wt - 0.25*(max(phi_vec) - min(phi_vec)), ...
				  phi_wt + 0.25*(max(phi_vec) - min(phi_vec)), 10), ...
		 ones(1,10)*targets(2), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	ylabel({'Net markup'; '(p90)'}, 'Interpreter', 'latex');
	xline(phi_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([phi_fix - .25, phi_fix + .25]);
	ylim([0,800])
	yticks([0,400,800])
	xticks([0.1,0.2,0.3,0.4,0.5])

	subplot(5,4,2)
	plot(lambda_vec(plot_vec), markup_p90_vec(2,plot_vec)-100, '--r', 'LineWidth', 2);
	hold on 
	plot(linspace(lambda - 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), ...
				  lambda + 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), 10), ...
		 ones(1,10)*targets(2), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	ylim([0, 200]);
	yticks([0,100,200]);
	xline(lambda_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([lambda_fix - .0155, lambda_fix + .0155]);
	xticks([1.01,1.02,1.03,1.04])


	subplot(5,4,3)
	plot((B_vec(plot_vec)*(1-0)^(-gamma))*12, markup_p90_vec(1,plot_vec)-100, ':k', 'LineWidth', 2);
	hold on 
	plot(linspace(12*(B*(1-0)^(-gamma)) - 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), ...
				  12*(B*(1-0)^(-gamma)) + 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), 10), ...
		 ones(1,10)*targets(2), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	ylim([0, 200]);
	yticks([0,100,200]);
	xline(B_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([0, 2.5]);
	xticks([0,.5,1,1.5,2,2.5])

	

	subplot(5,4,5)
	plot(phi_vec(plot_vec), RD_sales_quint_vec(3,plot_vec), '-b', 'LineWidth', 2);
	hold on 
	plot(linspace(phi_wt - 0.25*(max(phi_vec) - min(phi_vec)), ...
				  phi_wt + 0.25*(max(phi_vec) - min(phi_vec)), 10), ...
		 ones(1,10)*targets(3), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	ylabel({'R\&D to sales'; '(top profit quintile)'}, 'Interpreter', 'latex');
	xline(phi_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([phi_fix - .25, phi_fix + .25]);
	xticks([0.1,0.2,0.3,0.4,0.5])
	ylim([0,5])
	yticks([1,2,3,4,5])


	subplot(5,4,6)
	plot(lambda_vec(plot_vec),RD_sales_quint_vec(2,plot_vec),'--r', 'LineWidth', 2);
	hold on 
	plot(linspace(lambda - 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), ...
				  lambda + 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), 10), ...
		 ones(1,10)*targets(3), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(lambda_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([lambda_fix - .0155, lambda_fix + .0155]);
	ylim([0,5])
	yticks([1,2,3,4,5])
	xticks([1.01,1.02,1.03,1.04])

	subplot(5,4,7)
	plot((B_vec(plot_vec)*(1-0)^(-gamma))*12, RD_sales_quint_vec(1,plot_vec), ':k', 'LineWidth', 2);
	hold on 
	plot(linspace(12*(B*(1-0)^(-gamma)) - 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), ...
				  12*(B*(1-0)^(-gamma)) + 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), 10), ...
		 ones(1,10)*targets(3), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline((B_fix*(1-0)^(-gamma))*12, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([0,2.5]);
	ylim([0,5])
	yticks([1,2,3,4,5])
	xticks([0,.5,1,1.5,2,2.5])

	subplot(5,4,9)
	plot(phi_vec(plot_vec), innov_prof_vol_quint_vec(3,plot_vec), '-b', 'LineWidth', 2);
	hold on 
	plot(linspace(phi_wt - 0.25*(max(phi_vec) - min(phi_vec)), ...
				  phi_wt + 0.25*(max(phi_vec) - min(phi_vec)), 10), ...
		 ones(1,10)*targets(4), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$\phi$', 'Interpreter', 'latex');
	ylabel({'Profit vol.'; 'Top profit quintile'}, 'Interpreter', 'latex');
	xline(phi_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([phi_fix - .25, phi_fix + .25]);
	xticks([0.1,0.2,0.3,0.4,0.5])
	ylim([0,50])
	yticks([10,20,30,40,50])

	subplot(5,4,10)
	plot(lambda_vec(plot_vec),innov_prof_vol_quint_vec(2,plot_vec),'--r', 'LineWidth', 2);
	hold on 
	plot(linspace(lambda - 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), ...
				  lambda + 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), 10), ...
		 ones(1,10)*targets(4), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$\lambda$', 'Interpreter', 'latex');
	xline(lambda_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([lambda_fix - .0155, lambda_fix + .0155]);
	ylim([0,50])
	yticks([10,20,30,40,50])
	xticks([1.01,1.02,1.03,1.04])

	subplot(5,4,11)
	plot((B_vec(plot_vec)*(1-0)^(-gamma))*12, innov_prof_vol_quint_vec(1,plot_vec), ':k', 'LineWidth', 2);
	hold on 
	plot(linspace(12*(B*(1-0)^(-gamma)) - 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), ...
				  12*(B*(1-0)^(-gamma)) + 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), 10), ...
		 ones(1,10)*targets(4), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$B$', 'Interpreter', 'latex');
	xline((B_fix*(1-0)^(-gamma))*12, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([0, 2.5]);
	xticks([0,.5,1,1.5,2,2.5])
	ylim([0,50])
	yticks([10,20,30,40,50])

	subplot(5,4,13)
	plot(phi_vec(plot_vec), ppval(mod_store2{3,9}, phi_vec(plot_vec)), '-b', 'LineWidth', 2)
	hold on 
	plot(linspace(phi_wt - 0.25*(max(phi_vec) - min(phi_vec)), ...
				  phi_wt + 0.25*(max(phi_vec) - min(phi_vec)), 10), ...
		 ones(1,10)*targets(5), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$\phi$', 'Interpreter', 'latex');
	ylabel({'WITHIN'; ''}, 'Interpreter', 'latex');
	xline(phi_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([phi_fix - .25, phi_fix + .25]);
	xticks([0.1,0.2,0.3,0.4,0.5])
	ylim([75,100])
	yticks([75,80,85,90,95,100])

	subplot(5,4,14)
	plot(lambda_vec(plot_vec), ppval(mod_store2{2,9}, lambda_vec(plot_vec)), '--r', 'LineWidth', 2)
	hold on 
	plot(linspace(lambda - 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), ...
				  lambda + 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), 10), ...
		 ones(1,10)*targets(5), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$\lambda$', 'Interpreter', 'latex');
	xline(lambda_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([lambda_fix - .0155, lambda_fix + .0155]);
	ylim([75,100])
	yticks([75,80,85,90,95,100])
	xticks([1.01,1.02,1.03,1.04])

	subplot(5,4,15)
	plot((B_vec(plot_vec)*(1-0)^(-gamma))*12, ppval(mod_store2{1, 9}, B_vec(plot_vec)), ':k', 'LineWidth', 2)
	hold on 
	plot(linspace(12*(B*(1-0)^(-gamma)) - 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), ...
				  12*(B*(1-0)^(-gamma)) + 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), 10), ...
		 ones(1,10)*targets(5), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$B$', 'Interpreter', 'latex');
	xline((B_fix*(1-0)^(-gamma))*12, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([0,2.5]);
	xticks([0,.5,1,1.5,2,2.5])
	ylim([75,100])
	yticks([75,80,85,90,95,100])


	subplot(5,4,17)
	plot(phi_vec(plot_vec), ppval(mod_store2{3,8}, phi_vec(plot_vec)), '-b' ,'LineWidth', 2)
	hold on 
	plot(linspace(phi_wt - 0.25*(max(phi_vec) - min(phi_vec)), ...
				  phi_wt + 0.25*(max(phi_vec) - min(phi_vec)), 10), ...
		 ones(1,10)*targets(6), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$\phi$', 'Interpreter', 'latex');
	ylabel({'Innovation output'; '(p90)'}, 'Interpreter', 'latex');
	xline(phi_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([phi_fix - .25, phi_fix + .25]);
	xticks([0.1,0.2,0.3,0.4,0.5])
	ylim([0,40])
	yticks([0,10,20,30,40])
	xlabel('$\phi$', 'Interpreter', 'latex')

	subplot(5,4,18)
	plot(lambda_vec(plot_vec), ppval(mod_store2{2,8}, lambda_vec(plot_vec)), '--r', 'LineWidth', 2)
	hold on 
	plot(linspace(lambda - 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), ...
				  lambda + 0.25*(lambda_fix+0.015 - lambda_fix + 0.015), 10), ...
		 ones(1,10)*targets(6), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$\lambda$', 'Interpreter', 'latex');
	xline(lambda_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([lambda_fix - .0155, lambda_fix + .0155]);
	ylim([0,40])
	yticks([0,10,20,30,40])
	xticks([1.01,1.02,1.03,1.04])
	xlabel('$\lambda$', 'Interpreter', 'latex')

	subplot(5,4,19)
	plot((B_vec(plot_vec)*(1-0)^(-gamma))*12, ppval(mod_store2{1,8}, B_vec(plot_vec)), ':k', 'LineWidth', 2)
	hold on 
	plot(linspace(12*(B*(1-0)^(-gamma)) - 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), ...
				  12*(B*(1-0)^(-gamma)) + 0.25*((B_fix*(1-0)^(-gamma))*12 + 2 - ((B_fix*(1-0)^(-gamma))*12 - 2)), 10), ...
		 ones(1,10)*targets(6), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$B$', 'Interpreter', 'latex');
	xline((B_fix*(1-0)^(-gamma))*12, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([0, 2.5]);
	xticks([0,.5,1,1.5,2,2.5])
	ylim([0,40])
	yticks([0,10,20,30,40])
	xlabel('$B$', 'Interpreter', 'latex')


	subplot(5,4,4)
	plot(gamma_vec(plot_vec), markup_p90_vec(4,plot_vec)-100, '-b', 'color', [0, 0.5, 0], 'LineStyle', '-.', 'LineWidth', 2);
	hold on 
	plot(linspace(gamma - 0.25*(gamma_fix+0.18 - gamma_fix + 0.18), ...
				  gamma + 0.25*(gamma_fix+0.18 - gamma_fix + 0.18), 10), ...
		 ones(1,10)*targets(2), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(gamma_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([gamma_fix - .18, gamma_fix + .18]);
	xticks([0.3, 0.4,0.5,0.6])
	ylim([0,200])
	yticks([0,100,200])

	subplot(5,4,8)
	plot(gamma_vec(plot_vec), RD_sales_quint_vec(4,plot_vec), 'color', [0, 0.5, 0], 'LineStyle', '-.', 'LineWidth', 2);
	hold on 
	plot(linspace(gamma - 0.25*(gamma_fix+0.18 - gamma_fix + 0.18), ...
				  gamma + 0.25*(gamma_fix+0.18 - gamma_fix + 0.18), 10), ...
		 ones(1,10)*targets(3), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xline(gamma_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([gamma_fix - .18, gamma_fix + .18]);
	xticks([0.3,0.4,0.5,0.6])
	ylim([0,5])
	yticks([1,2,3,4,5])

	subplot(5,4,12)
	plot(gamma_vec(plot_vec), innov_prof_vol_quint_vec(4,plot_vec), 'color', [0, 0.5, 0], 'LineStyle', '-.', 'LineWidth', 2);
	hold on 
	plot(linspace(gamma - 0.25*(gamma_fix+0.18 - gamma_fix + 0.18), ...
				  gamma + 0.25*(gamma_fix+0.18 - gamma_fix + 0.18), 10), ...
		 ones(1,10)*targets(4), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$\gamma$', 'Interpreter', 'latex');
	xline(gamma_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([gamma_fix - .18, gamma_fix + .18]);
	xticks([0.3,0.4,0.5,0.6])
	ylim([0,50])
	yticks([10,20,30,40,50])


	subplot(5,4,16)
	plot(gamma_vec(plot_vec), ppval(mod_store{4,9}, gamma_vec(plot_vec)), 'color', [0,0.5,0], 'LineStyle', '-.', 'LineWidth', 2)
	hold on 
	plot(linspace(gamma - 0.25*(gamma_fix+0.18 - gamma_fix + 0.18), ...
				  gamma + 0.25*(gamma_fix+0.18 - gamma_fix + 0.18), 10), ...
		 ones(1,10)*targets(5), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$\gamma$', 'Interpreter', 'latex');
	xline(gamma_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([gamma_fix - .18, gamma_fix + .18]);
	xticks([0.3,0.4,0.5,0.6])
	ylim([75,100])
	yticks([75,80,85,90,95,100])

	subplot(5,4,20)
	plot(gamma_vec(plot_vec), ppval(mod_store2{4,8}, gamma_vec(plot_vec)), 'color', [0, 0.5, 0], 'LineStyle', '-.', 'LineWidth', 2)
	hold on 
	plot(linspace(gamma - 0.25*(gamma_fix+0.18 - gamma_fix + 0.18), ...
				  gamma + 0.25*(gamma_fix+0.18 - gamma_fix + 0.18), 10), ...
		 ones(1,10)*targets(6), ...
		 'LineStyle', ':', 'color', [0, 0.5, 0]);
	xlabel('$\gamma$', 'Interpreter', 'latex');
	xline(gamma_fix, 'LineStyle', '-', 'color', [.5,.5,.5], 'LineWidth', 1.25);
	xlim([gamma_fix - .18, gamma_fix + .18]);
	xticks([0.3, 0.4,0.5,0.6])
	ylim([0,40])
	yticks([0,10,20,30,40])
	xlabel('$\gamma$', 'Interpreter', 'latex')


	saveas(gcf, strcat("Output/Figures_Paper/fourmom_text_IA.eps"), 'epsc');
	saveas(gcf, strcat("Output/Figures_Paper/fourmom_text_IA.png"));








end