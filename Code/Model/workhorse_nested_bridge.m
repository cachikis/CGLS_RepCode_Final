%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: workhorse_nested_bridge.m
% Author: Craig A. Chikis
% Date: 12/14/2020
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function return_val = workhorse_nested_bridge(s_max, w_guess, v_guess, rho, xbar, eta, lab_tol, val_tol, ...
		  										 B, gamma, lambda, entrance, subsidy, tax_rate, ...
		  										 subsidyE, B_entrant, PE, entry_type, ...
		  										 prob, probE, prob_exog, ...
		  										 alpha_fric, profit, EIS, r_guess, find_r, ...
		  										 elastic_labor, v_max_try, nu_s, kappa)


probF = flip(prob(1:s_max, :), 1);
probL = prob((s_max+2):end, :);
prob0 = prob(s_max+1, :);
probE = flip(probE, 1);
prob_exog = flip(prob_exog, 1);

eta(1) = 0;

value_functions = nan(1, 2*s_max+1);
investment = nan(1, 2*s_max+1);
firm_distributions = nan(1, s_max+1);
growth_rate = nan;
markup = nan;
wage_share = nan;
RMSE = nan;
value_flag = false;
labor_flag = false;
welfare = nan;
investment_entrant = nan(1, s_max+1);

clear_labor = false;
clear_value = false;

count_labor = 0;
max_labor = 1e10;
max_value = v_max_try;
l_v_condition = 1e10; 

labor_diff = 1e10;
val_diff = 1e10;

G = @(x) (x/B).^(1/gamma);
G_entrant = @(x) (x/B_entrant).^(1/gamma);
GPInv = @(z) B*((z*B*gamma).^(gamma/(1-gamma)));
GPInv_entrant = @(z) B_entrant*((z*B_entrant*gamma).^(gamma/(1-gamma)));
GInv = @(z,B) B*z.^gamma;


xE = zeros(1, s_max+1);

mu_temp = repelem(1/length(0:s_max), length(0:s_max));
ub_w = 1;
lb_w = 0;

	count_value = 0;
	clear_value = false;
	value_flag = false;
	l_v_condition = 1e10;


	if (elastic_labor)
		w = 1;
	else
		w = r_guess(2);
	end

	while (~clear_value)
		count_value = count_value + 1;

		vL = v_guess((s_max+2):end);
		vF = flip(v_guess(1:s_max));
		v0 = v_guess(s_max+1);

		FOC_L = zeros(1, s_max);
		FOC_F = zeros(1, s_max);

		for (j = 1:length(1:s_max))
			probL_it = probL(j, :);
			FOC_L(j) = sum(probL_it.*v_guess) - vL(j);

			probF_it = probF(j, :);
			FOC_F(j) = sum(probF_it.*v_guess) - vF(j);

		end

		FOC_0 = sum(prob0.*v_guess) - v0;


		xL = max(min(max(GPInv(FOC_L./((1-subsidy((s_max+2):end))*w)), 0), ...
					 GInv(alpha_fric*vL./(w*(1-subsidy((s_max+2):end))), B)), 0);

		xF = max(min(max(GPInv(FOC_F./((1-flip(subsidy(1:s_max)))*w)), 0), ...
					 GInv(alpha_fric*vF./(w*(1-flip(subsidy(1:s_max)))), B)), 0);

		x0 = max(min(max(GPInv(FOC_0/((1-subsidy(s_max+1))*w)), 0), GInv(alpha_fric*v0/(w*(1-subsidy(s_max+1))), B)),0);

		if (entrance)

			if (strcmp(entry_type, "Ufuk/Sina")) 
			
				FOC_E = zeros(1, length(0:s_max));
				probE_tmp = probE;
				for (j = 1:length(0:s_max))
					FOC_E(j) = sum(probE_tmp(j,:).*v_guess) - 0;
				end

				xE = max(GPInv_entrant(FOC_E./((1-subsidyE)*w)), 0);

			elseif (strcmp(entry_type, "Undirected"))
				if (mod(count_value, 50) == 0) && ((count_labor > 1) || (count_value > 1)) 

					inv_temp = [flip(xF), x0, xL];
					inv_temp_E = xE;

					[A, mu_temp, b] = firm_distributions_v3(inv_temp, inv_temp_E, eta, probL, probF, prob0, ...
												 probE, prob_exog, s_max);

					mu_temp = mu_temp';

				end

				
				
				FOC_E = zeros(1, s_max+1);
				probE_tmp = probE;
				for (j = 1:length(0:s_max))
					FOC_E(j) = sum(probE_tmp(j,:).*v_guess) - 0;
				end

				FOC_E = sum(mu_temp.*FOC_E);
				wage_div = (1/w)*(1/(1-subsidyE(2)));

				
				xE = max(GPInv_entrant(FOC_E*wage_div), 0);
				xE = repelem(xE, length(0:s_max));


			end
			
		end


		xL(imag(xL) ~= 0) = 0;
		xF(imag(xF) ~= 0) = 0;
		x0(imag(x0) ~= 0) = 0;
		xE(imag(xE) ~= 0) = 0;

		eta_max = max(eta);

		beta = (2*xbar + eta_max)/(rho + 2*xbar + eta_max)*(~entrance) + ...
			   (3*xbar + eta_max)/(rho + 3*xbar + eta_max)*entrance; 

		beta = repelem(beta, length(0:s_max));

		value_function_in = [flip(vF), v0, vL];
		investment_in = [flip(xF), x0, xL];



		wage_share = w;
		v_guess = vfi_main(value_function_in, investment_in, eta, ...
													beta, s_max, ...
													wage_share, B, gamma, xbar, entrance, ...
													subsidy, lambda, rho, xE, ...
													B_entrant, probL, probF, prob0, ...
													probE, prob_exog, profit, EIS, r_guess(1));


		
		v_condition = norm(v_guess - value_function_in)/sum(v_guess);




		l_v_condition = v_condition;
		RMSE = mean(norm(v_guess - value_function_in))^.5;

		if (v_condition < val_tol)
			clear_value = true;
			value_flag = true;
		elseif (count_value >= max_value)
			clear_value = true;
		end

		if (xbar < max(investment_in))
			xbar = max(investment_in) + 5;
		end


	end


	[A, mu, b] = firm_distributions_v3(investment_in, xE, eta, probL, probF, prob0, ...
								    probE, prob_exog, s_max);
	mu = mu';

	L_demand = [G(x0), G(xL)];
	F_demand = [G(x0), G(xF)];
	E_demand = G_entrant(xE);
	if (kappa >= 9999) 
		non_RD_demand = lambda.^(-(0:s_max)) * (1/w);
	else
		non_RD_demand = (nu_s.^(1-kappa))*(kappa-1)*(1/w)./((1+nu_s.^(1-kappa)).*(kappa+nu_s.^(1-kappa))) + ...
						(nu_s.^(1-kappa))*(kappa-1)*(1/w)./((1+nu_s.^(1-kappa)).*(kappa*nu_s.^(1-kappa)+1));
	end

	labor_diff = 1 - sum(mu.*(L_demand + F_demand + E_demand + non_RD_demand));

		clear_labor = true;
		labor_flag = true;
		value_functions = value_function_in;
		investment = [flip(xF), x0, xL];
		firm_distributions = mu;
		if (kappa >= 9999)
			markup = sum(mu.*(lambda.^(0:s_max) - 1)); 
		else
			markup = nan;
		end
		wage_share = w;
		investment_entrant = xE;



		gE_iter = zeros(1, length(0:s_max));
		gL_iter = zeros(1, length(0:s_max));
		gF_iter = zeros(1, length(0:s_max));

		probL_iter = [prob0;probL];
		probF_iter = [prob0;probF];
		

		for (s = 0:s_max)
			gE_iter(s+1) = sum(probE(s+1,(s_max+2):end).*(1:s_max));
			gL_iter(s+1) = sum(probL_iter(s+1,(s_max+2+s):end).*(((1+s):s_max)-s));
			gF_iter(s+1) = sum(probF_iter(s+1, (s_max+2):end).*(1:s_max));
		end

		growth_rate = log(lambda)*sum(mu.*(xE.*gE_iter + [x0, xF].*gF_iter + [x0, xL].*gL_iter));


		return_val = {value_functions, investment, firm_distributions, ...
					  growth_rate, markup, wage_share, RMSE, ...
					  value_flag, labor_flag, ...
					  investment_entrant};


if (find_r)
	return_val = zeros(size(r_guess));
	if (elastic_labor)
		return_val = growth_rate - EIS*(r_guess - rho);
	else
		return_val(1) = growth_rate - EIS*(r_guess(1) - rho);
		return_val(2) = labor_diff;
	end
end


end
