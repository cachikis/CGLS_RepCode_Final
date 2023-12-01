%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: find_guess_fsolve.m
% Author: Craig A. Chikis
% Date: 11/30/2020
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [v_guess, err, err_rho, wsolve] = find_guess_fsolve(v_guess, prob, subsidy, w_guess, elastic_labor, profit, ...
							    	 probE, prob_exog, B, B_entrant, lambda, gamma, subsidyE, ...
							         tax_rate, xbar, entrance, entry_type, eta, rho, val_tol, ...
							         alpha_fric, EIS, r_guess, v_max_try, LMS_bar_s)



function out = search_w(inp, return_criterion, v_guess, prob, subsidy, elastic_labor, profit, ...
						probE, prob_exog, B, B_entrant, lambda, gamma, subsidyE, ...
						tax_rate, xbar, entrance, entry_type, eta, rho, val_tol, ...
						alpha_fric, EIS, find_r, v_max_try, LMS_bar_s, r_guess)

	w = inp(1);
	try 
		r_guess = inp(2);
	catch
	end 

	count_value = 0;
	s_max = (length(v_guess)-1)/2;
	clear_value = false;
	max_value = v_max_try;

	G = @(x,B) (x/B).^(1/gamma);
	GPInv = @(z,B) B*((z*B*gamma).^(gamma/(1-gamma)));
	GInv = @(z,B) B*z.^gamma;

	probL = prob((s_max+2):end, :);
	probF = flip(prob(1:s_max, :), 1);
	prob0 = prob(s_max+1, :);
	probE = flip(probE,1);
	prob_exog = flip(prob_exog,1);




	mu_temp = repelem(1/length(0:s_max), length(0:s_max));
	xE = zeros(1, s_max+1);
	while (~clear_value)
		count_value = count_value + 1;

		vL = v_guess((s_max+2):end);
		vF = flip(v_guess(1:s_max));
		v0 = v_guess(s_max+1);

		FOC_L = zeros(1, s_max);
		FOC_F = zeros(1, s_max);

		for (ii = 1:length(1:s_max))
			probL_it = prob(s_max+1+ii, :);
			FOC_L(ii) = sum(probL_it.*v_guess) - v_guess(s_max+1+ii);

			probF_it = prob(s_max+1-ii, :);
			FOC_F(ii) = sum(probF_it.*v_guess) - v_guess(s_max+1-ii);
		end
		FOC_0 = sum(prob(s_max+1,:).*v_guess) - v_guess(s_max+1);

		xL = max(GPInv(FOC_L./((1-subsidy((s_max+2):end))*w), B), 0);
		xF = max(GPInv(FOC_F./((1-flip(subsidy(1:s_max)))*w), B), 0);
		x0 = max(GPInv(FOC_0/((1-subsidy(s_max+1))*w), B), 0);

		if (entrance)

			if (strcmp(entry_type, "Ufuk/Sina"))
					
				FOC_E = zeros(1, length(0:s_max));
				for (ii = 1:length(0:s_max))
					FOC_E(ii) = sum(probE(ii,:).*v_guess) - 0;
				end

				xE = max(GPInv(FOC_E./((1-subsidyE)*w), B_entrant), 0);

			elseif (strcmp(entry_type, "Undirected"))
				if (mod(count_value, 50) == 0) && (count_value > 1)

					inv_temp = [flip(xF), x0, xL];
					inv_temp_E = xE;

					[A, mu_temp, b] = firm_distributions_v3(inv_temp, inv_temp_E, eta, probL, probF, prob0, ...
															probE, prob_exog, s_max);

					mu_temp = mu_temp';

				end
				FOC_E = zeros(1, s_max+1);
				for (ii = 1:length(0:s_max))
					FOC_E(ii) = sum(probE(ii,:).*v_guess) - 0;
				end

				FOC_E = sum(mu_temp.*FOC_E);
				wage_div = (1/w)*(1/(1-subsidyE(2)));

			
				xE = max(GPInv(FOC_E*wage_div, B_entrant), 0);
				xE = repelem(xE, length(0:s_max));


			end
					
		end

		xL(imag(xL) ~= 0) = 0;
		xF(imag(xF) ~= 0) = 0;
		x0(imag(x0) ~= 0) = 0;
		xE(imag(xE) ~= 0) = 0;

		eta_max = max(eta);
		rho_term = EIS*rho + (1-EIS)*r_guess;


		beta = ((2*xbar + eta)./(rho_term + 2*xbar + eta)).*(~entrance) + ...
			((3*xbar + eta)./(rho_term + 3*xbar + eta)).*entrance; 

		

		value_function_in = [flip(vF), v0, vL];
		investment_in = [flip(xF), x0, xL];

		investment_in = max(min(investment_in, GInv(alpha_fric*value_function_in./(w*(1-subsidy)), B)), 0);
		investment_in(imag(investment_in) ~= 0) = 0;

		wage_share = w;
		v_guess = vfi_main(value_function_in, investment_in, eta, ...
														beta, s_max, ...
														wage_share, B, gamma, xbar, entrance, ...
														subsidy, lambda, rho, xE, ...
														B_entrant, probL, probF, prob0, ...
														probE, prob_exog, profit, EIS, r_guess);

				
		v_condition = norm(v_guess - value_function_in)/sum(v_guess);


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
	err = sqrt(mean((v_guess - value_function_in).^2));
	err_rho = sqrt(mean((rho*v_guess - rho*value_function_in).^2));

	[A, mu, b] = firm_distributions_v3(investment_in, xE, eta, probL, probF, prob0, ...
															probE, prob_exog, s_max);
	mu = mu';

	if (elastic_labor) && (entrance)
		v_guess = [v_guess, mu];
	elseif (~elastic_labor) && (entrance)
		v_guess = [v_guess, mu, w];
	elseif (elastic_labor) && (~entrance)
		v_guess = v_guess;
	elseif (~elastic_labor) && (~entrance)
		v_guess = [v_guess, w];
	end

	if (EIS ~= 1) 
		v_guess = [v_guess, r_guess];
	end


	if (EIS ~= 1) && (find_r)
		gE_iter = zeros(1, length(0:s_max));
		gL_iter = zeros(1, length(0:s_max));
		gF_iter = zeros(1, length(0:s_max));

		probL_iter = prob((s_max+1):end,:);
		probF_iter = flip(prob(1:(s_max+1), :),1);
				

		for (s = 0:s_max)
			gE_iter(s+1) = sum(probE(s+1,(s_max+2):end).*(1:s_max));
			gL_iter(s+1) = sum(probL_iter(s+1,(s_max+2+s):end).*(((1+s):s_max)-s));
			gF_iter(s+1) = sum(probF_iter(s+1, (s_max+2):end).*(1:s_max));
		end

		growth_rate = log(lambda)*sum(mu.*(xE.*gE_iter + ...
										flip(investment_in(1:(s_max+1))).*gF_iter + ...
										investment_in((s_max+1):end).*gL_iter));

		v_guess = growth_rate - EIS*(r_guess - rho);
	end

	lab_diff = 1 - sum(mu.*(G([x0, xL], B) + G([x0, xF], B) + ...
								G(xE, B_entrant) + 1./(w*min(LMS_bar_s, lambda.^(0:s_max)))));

	if (~return_criterion)
		out = {v_guess, err, err_rho}; 
	else
		if (EIS ~= 1)
			out = [lab_diff, v_guess]';
		else
			out = lab_diff; 
		end
	end

end


if (elastic_labor) && (EIS == 1)
	return_criterion = false;
	inp = [1, r_guess]'; 
	find_r = false;
	out = search_w(inp, return_criterion, v_guess, prob, subsidy, elastic_labor, profit, ...
						probE, prob_exog, B, B_entrant, lambda, gamma, subsidyE, ...
						tax_rate, xbar, entrance, entry_type, eta, rho, val_tol, ...
						alpha_fric, EIS, find_r, v_max_try, LMS_bar_s, r_guess);
	v_guess = out{1}; 
	err = out{2}; 
	err_rho = out{3};
	wsolve = inp(1);
else
	return_criterion = true;
	find_r = true; 
	options = optimset('Display', 'off', 'MaxIter', 1e4, 'MaxFunEval', 1e4, ...
				   	   'TolX', eps, 'TolFun', eps, 'Algorithm', 'trust-region-dogleg', ...
				        'UseParallel', false);	
	options_2 = optimset('Display', 'off'); 
	inp = [w_guess, r_guess]'; 
	outfunc = @(inp) search_w(inp, return_criterion, v_guess, prob, subsidy, elastic_labor, profit, ...
						probE, prob_exog, B, B_entrant, lambda, gamma, subsidyE, ...
						tax_rate, xbar, entrance, entry_type, eta, rho, val_tol, ...
						alpha_fric, EIS, find_r, v_max_try, LMS_bar_s, r_guess);
	if (EIS ~= 1)
		wsolve = fsolve(outfunc, inp, options); 
	else
		wsolve = fzero(outfunc, inp(1), options_2); 
	end
	return_criterion = false; 
	find_r = false; 
	out = search_w(wsolve, return_criterion, v_guess, prob, subsidy, elastic_labor, profit, ...
						probE, prob_exog, B, B_entrant, lambda, gamma, subsidyE, ...
						tax_rate, xbar, entrance, entry_type, eta, rho, val_tol, ...
						alpha_fric, EIS, find_r, v_max_try, LMS_bar_s, r_guess);
	v_guess = out{1}; 
	err = out{2}; 
	err_rho = out{3};
	wsolve = wsolve(1);
end

end


