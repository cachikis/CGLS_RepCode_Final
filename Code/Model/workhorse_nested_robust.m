%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: workhorse_nested_robust.m
% Author: Craig A. Chikis
% Date: 11/30/2020
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = workhorse_nested_robust(m)



	s_max = m.s_max;
	w_guess = m.w_guess;
	v_guess = m.v_guess;
	rho = m.rho;
	xbar = m.xbar;
	eta = m.eta;
	lab_tol = m.lab_tol;
	val_tol = m.val_tol; 
	B = m.B;
	gamma = m.gamma; 
	lambda = m.lambda; 
	entrance = m.entrance; 
	subsidy = m.subsidy; 
	tax_rate = m.tax_rate;
	subsidyE = m.subsidyE;
	B_entrant = m.B_entrant;
	PE = m.PE;
	entry_type = m.entry_type;
	prob = m.prob; 
	probE = m.probE; 
	prob_exog = m.prob_exog; 
	kappa = m.kappa; 
	LMS_bar_s = m.LMS_bar_s; 
	elastic_labor = m.elastic_labor; 
	algorithm_fsolve = m.algorithm_fsolve; 
	alpha_fric = m.alpha_fric; 
	bridge_method = m.bridge_method; 
	EIS = m.EIS; 
	r_guess = m.r_guess; 
	v_max_try = m.v_max_try; 
	two_adjust = m.two_adjust; 

	if (EIS ~= 1) && (elastic_labor)
		error('This code cannot run with nonunity EIS.')
	elseif (EIS ~= 1) && (~elastic_labor) && (~(kappa >= 9999))
		error('Not sure what to do with w for calculating production labor for labor market clearing')
	end
	
	
	[profit, nu_s] = get_profit(s_max, lambda, kappa, tax_rate, LMS_bar_s, two_adjust);

	options = optimset('Display', 'off', 'MaxIter', 5e4, 'MaxFunEval', 5e4, ...
					'TolX', eps, 'TolFun', eps, 'Algorithm', algorithm_fsolve, ...
					'UseParallel', false);

	options_r = optimset('Display', 'off');



	if (all(prob(1:(s_max-1), s_max+1) == 0)) && (~elastic_labor) && (kappa >= 9999) && ((abs(eta(2) - 0.024/12) < 1e-8) || ...
		(abs(eta(2) - 0.5*0.024/12) < 1e-8) || (abs(eta(2) - 0.0005) < 1e-8))   && ...
			(EIS == 1)
		[value_functions, investment, firm_distributions, ...
				growth_rate, markup, wage_share, RMSE, ...
				value_flag, labor_flag, ...
				investment_entrant] = workhorse_nested(s_max, w_guess, v_guess, rho, xbar, eta, lab_tol, val_tol, ...
														B, gamma, lambda, entrance, subsidy, tax_rate, ...
														subsidyE, B_entrant, PE, entry_type, ...
														prob, probE, prob_exog, v_max_try);
	else
		exitflag_r = 1;

		[v_in0, vfi_error, vfi_error_rho, wsolve] = find_guess_fsolve(v_guess, prob, subsidy, w_guess, elastic_labor, profit, ...
										probE, prob_exog, B, B_entrant, lambda, gamma, subsidyE, ...
										tax_rate, xbar, entrance, entry_type, eta, rho, val_tol, ...
										alpha_fric, EIS, r_guess, v_max_try, LMS_bar_s);

				
		solve_Upsilon = @(v_in) satisfy_Upsilon(rho, v_in, prob, probE, prob_exog, subsidy, subsidyE, ...
											tax_rate, B, gamma, entrance, B_entrant, lambda, eta, ...
											profit, elastic_labor, entry_type, s_max, ...
											alpha_fric, EIS, kappa, LMS_bar_s);

		vfi_hjb_error = sqrt(mean(solve_Upsilon(v_in0).^2));
		vfi_hjb_error_rho = sqrt(mean((rho*solve_Upsilon(v_in0)).^2));



		[v_in1, fval1, flag1, output] = fsolve(solve_Upsilon, v_in0, options);
		flag2 = 0;
		end_hjb_error = sqrt(mean(solve_Upsilon(v_in1).^2));
		end_hjb_error_rho = sqrt(mean((rho*solve_Upsilon(v_in1)).^2));
			
	
		if (~(flag1 > 0)) && (~strcmp(algorithm_fsolve, 'levenberg-marquardt')) 
			options = optimset('Display', 'off', 'MaxIter', 5e4, 'MaxFunEval', 5e4, ...
							'TolX', eps, 'TolFun', eps, 'Algorithm', 'levenberg-marquardt', ...
							'UseParallel', false);

			
			[v_in2, fval2, flag2, output] = fsolve(solve_Upsilon, v_in0, options);
			end_hjb_error = sqrt(mean(solve_Upsilon(v_in2).^2));
			end_hjb_error_rho = sqrt(mean((rho*solve_Upsilon(v_in2)).^2));


		end
		flag = max(flag1, flag2);


		if (((all(prob(1:(s_max-1), s_max+1) == 0)) && ...
													(~elastic_labor) && (~(flag > 0)) && (EIS == 1) ...
														)) && ...
					(bridge_method) 
			flag = 1;
			if (sum(abs(fval1)) < sum(abs(fval2)))
				v_in = v_in1;
				fval = fval1;
			else
				v_in = v_in2;
				fval = fval2;
			end

			if (~elastic_labor) && ((v_in(end) > 1 || v_in(end) < 0))
				flag = 0;
			end

		else
			if (flag1 > 0)
				v_in = v_in1;
				fval = fval1;
			elseif (flag2 > 0)
				v_in = v_in2;
				fval = fval2;
			else
				v_in = v_in1;
				fval = fval1;
			end
		end


		exitflag_r2 = 1;
		if (~(flag > 0)) && (bridge_method) 
			if (elastic_labor)
				w_guess = 1;
				lab_tol = Inf;
			end
			v_guess_fallback = v_in(1:(2*s_max+1));

				if (elastic_labor)
					find_r2 = true;
					r_guess2 = r_guess;
					unknown_r2 = @(r_guess) workhorse_nested_bridge(s_max, w_guess, v_guess, rho, xbar, eta, lab_tol, val_tol, ...
														B, gamma, lambda, entrance, subsidy, tax_rate, ...
														subsidyE, B_entrant, PE, entry_type, ...
														prob, probE, prob_exog, ...
														alpha_fric, profit, EIS, r_guess, find_r2, elastic_labor, ...
														v_max_try, nu_s, kappa);

					[found_r2, ~, exitflag_r2] = fzero(unknown_r2, r_guess2, options_r);
				else
					find_r2 = true;
					r_guess2 = [r_guess, w_guess];

					unknown_r2 = @(r_guess) workhorse_nested_bridge(s_max, w_guess, v_guess_fallback, ...
														rho, xbar, eta, lab_tol, val_tol, ...
														B, gamma, lambda, entrance, subsidy, tax_rate, ...
														subsidyE, B_entrant, PE, entry_type, ...
														prob, probE, prob_exog, ...
														alpha_fric, profit, EIS, r_guess, find_r2, elastic_labor, ...
														v_max_try, nu_s, kappa);
					options = optimset('Display', 'off', 'MaxIter', 50000, 'MaxFunEval', 50000, ...
									'TolX', eps, 'TolFun', eps, 'Algorithm', algorithm_fsolve, ...
									'UseParallel', false);
					[found_r2, ~, exitflag_r2] = fsolve(unknown_r2, r_guess2, options);
				end
				find_r = false;
				return_bridge = workhorse_nested_bridge(s_max, w_guess, v_guess, rho, xbar, eta, lab_tol, val_tol, ...
													B, gamma, lambda, entrance, subsidy, tax_rate, ...
													subsidyE, B_entrant, PE, entry_type, ...
													prob, probE, prob_exog, ...
													alpha_fric, profit, EIS, found_r2, find_r, elastic_labor, ...
													v_max_try, nu_s, kappa);

			

			value_functions_t = return_bridge{1};
			investment_t = return_bridge{2};
			firm_distributions_t = return_bridge{3};
			growth_rate_t = return_bridge{4};
			markup_t = return_bridge{5};
			wage_share_t = return_bridge{6};
			RMSE_t = return_bridge{7};
			value_flag_t = return_bridge{8};
			labor_flag_t = return_bridge{9};
			investment_entrant_t = return_bridge{10};

			if (elastic_labor) && (EIS == 1) && (entrance)
				v_in = [value_functions_t, firm_distributions_t];
			elseif (elastic_labor) && (EIS ~= 1) && (entrance)
				v_in = [value_functions_t, firm_distributions_t, found_r2];
			elseif (~elastic_labor) && (EIS == 1) && (entrance)
				v_in = [value_functions_t, firm_distributions_t, wage_share_t];
			elseif (~elastic_labor) && (EIS ~= 1) && (entrance)
				v_in = [value_functions_t, firm_distributions_t, wage_share_t, found_r2(1)];
			end

			if (elastic_labor) && (EIS == 1) && (~entrance)
				v_in = [value_functions_t];
			elseif (elastic_labor) && (EIS ~= 1) && (~entrance)
				v_in = [value_functions_t, found_r2];
			elseif (~elastic_labor) && (EIS == 1) && (~entrance)
				v_in = [value_functions_t, wage_share_t];
			elseif (~elastic_labor) && (EIS ~= 1) && (~entrance)
				v_in = [value_functions_t, wage_share_t, found_r2(1)];
			end



			flag = 1;

		end
	

		if (flag > 0) && (exitflag_r > 0) && (exitflag_r2 > 0) 


			[value_functions, investment, firm_distributions, ...
						growth_rate, markup, wage_share, RMSE, ...
						value_flag, labor_flag, ...
						investment_entrant, ...
						markup_eq, markup_eq_array, ...
						markup_wt, markup_wt_array, ...
						lab_demand, ...
						real_interest_rate, ...
						growth_rate_LMS, ...
						fric_binding] = get_BGP(v_in, elastic_labor, prob, probE, prob_exog, ...
																subsidy, lambda, nu_s, rho, B, B_entrant, gamma, ...
																entrance, entry_type, subsidyE, kappa, alpha_fric, ...
																EIS, LMS_bar_s, eta);


					
		end
	end
	try
		m.value_functions = value_functions;
	catch
		m.value_functions = nan(1,2*s_max+1);
	end
	try 
		m.investment = investment;
	catch
		m.investment = nan(1,2*s_max+1);
	end
	try
		m.firm_distributions = firm_distributions;
	catch
		m.firm_distributions = nan(1,s_max+1); 
	end
	try
		m.growth_rate = growth_rate; 
	catch
		m.growth_rate = nan;
	end
	try
		m.markup = markup; 
	catch
		m.markup = nan;
	end
	try
		m.wage_share = wage_share;
	catch
		m.wage_share = nan;
	end
	try
		m.RMSE = RMSE; 
	catch
		m.RMSE = nan;
	end
	try
		m.value_flag = value_flag; 
	catch
		m.value_flag = nan;
	end
	try
		m.labor_flag = labor_flag; 
	catch
		m.labor_flag = nan;
	end
	try 
		m.investment_entrant = investment_entrant;
	catch
		m.investment_entrant = nan(1,s_max+1);
	end
	try
		m.markup_eq = markup_eq; 
	catch
		m.markup_eq = nan;
	end
	try
		m.markup_eq_array = markup_eq_array;
	catch
		m.markup_eq_array = nan(1,s_max+1);
	end
	try 
		m.markup_wt = markup_wt;
	catch
		m.markup_wt = nan;
	end
	try
		m.markup_wt_array = markup_wt_array;
	catch
		m.markup_wt_array = nan(1,s_max+1);
	end
	try
		m.profit = profit;
	catch
		m.profit = nan(1,2*s_max+1);
	end
	try
		m.nu_s = nu_s; 
	catch
		m.nu_s = nan(1,2*s_max+1); 
	end
	try
		m.lab_demand = lab_demand; 
	catch
		m.lab_demand = cell(1,3);
	end
	try
		m.real_interest_rate = real_interest_rate;
	catch
		m.real_interest_rate = nan;
	end
	try
		m.vfi_error = vfi_error;
	catch
		m.vfi_error = nan;
	end
	try
		m.vfi_hjb_error = vfi_hjb_error;
	catch
		m.vfi_hjb_error = nan;
	end
	try
		m.end_hjb_error = end_hjb_error;
	catch
		m.end_hjb_error = nan;
	end
	try
		m.vfi_error_rho = vfi_error_rho;
	catch
		m.vfi_error_rho = nan;
	end
	try
		m.vfi_hjb_error_rho = vfi_hjb_error_rho;
	catch
		m.vfi_hjb_error_rho = nan;
	end
	try
		m.end_hjb_error_rho = end_hjb_error_rho;
	catch
		m.end_hjb_error_rho = nan;
	end
	try 
		m.fric_binding = fric_binding;
	catch
		m.fric_binding = nan(1, 2*m.s_max+1); 
	end



end





