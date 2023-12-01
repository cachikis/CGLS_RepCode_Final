%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: satisfy_Upsilon.m
% Author: Craig A. Chikis
% Note(s):
% MA-MFA
% Date: 11/24/2020
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Upsilon = satisfy_Upsilon(rho, v_in, prob, probE, prob_exog, subsidy, subsidyE, ...
								   tax_rate, B, gamma, entrance, B_entrant, lambda, eta, ...
								   profit, elastic_labor, entry_type, s_max, alpha_fric, ...
								   EIS, kappa, LMS_bar_s)
	
	G = @(x,B) (x/B).^(1/gamma);
	GPInv = @(z,B) B*((z*B*gamma).^(gamma/(1-gamma)));
	GInv = @(z,B) B*z.^gamma;

	probL = prob((s_max+2):end, :);
	probF = flip(prob(1:s_max, :), 1);
	prob0 = prob(s_max+1, :);
	probE = flip(probE,1);
	prob_exog = flip(prob_exog,1);

	v_guess = v_in(1:(2*s_max+1));
	if (entrance)
		mu_guess = v_in((2*s_max+2):(3*s_max+2));
	end

	if (elastic_labor) && (EIS == 1) && (entrance)
		Upsilon = nan(1, 3*s_max+2);
		w_guess = 1; 
		r_guess = 0;
	elseif (~elastic_labor) && (EIS == 1) && (entrance)
		Upsilon = nan(1, 3*s_max+3);
		w_guess = v_in(3*s_max+3);
		r_guess = 0;
	elseif (elastic_labor) && (EIS ~= 1) && (entrance)
		Upsilon = nan(1, 3*s_max+3);
		w_guess = 1;
		r_guess = v_in(3*s_max+3);
	elseif (~elastic_labor) && (EIS ~= 1) && (entrance)
		Upsilon = nan(1, 3*s_max+4);
		w_guess = v_in(3*s_max+3);
		r_guess = v_in(3*s_max+4);
	end

	if (elastic_labor) && (EIS == 1) && (~entrance)
		Upsilon = nan(1, 2*s_max+1);
		w_guess = 1; 
		r_guess = 0;
	elseif (~elastic_labor) && (EIS == 1) && (~entrance)
		Upsilon = nan(1, 2*s_max+2);
		w_guess = v_in(2*s_max+2);
		r_guess = 0;
	elseif (elastic_labor) && (EIS ~= 1) && (~entrance)
		Upsilon = nan(1, 2*s_max+2);
		w_guess = 1;
		r_guess = v_in(2*s_max+2);
	elseif (~elastic_labor) && (EIS ~= 1) && (~entrance)
		Upsilon = nan(1, 2*s_max+3);
		w_guess = v_in(2*s_max+2);
		r_guess = v_in(2*s_max+3);
	end

	FOC = zeros(1, 2*s_max+1);
	for (ii = 1:s_max)
		FOC(s_max+1+ii) = sum(prob(s_max+1+ii, :).*v_guess) - v_guess(s_max+1+ii);
		FOC(s_max+1-ii) = sum(prob(s_max+1-ii, :).*v_guess) - v_guess(s_max+1-ii);
	end
	FOC(s_max+1) = sum(prob(s_max+1, :).*v_guess) - v_guess(s_max+1);

	investment_in = max(GPInv(FOC./(w_guess*(1-subsidy)), B), 0);
	investment_in(imag(investment_in) ~= 0) = 0;	

	investment_in = max(min(investment_in, GInv(alpha_fric.*v_guess./(w_guess*(1-subsidy)), B)), 0);
	investment_in(imag(investment_in) ~= 0) = 0;

	xE = zeros(1, s_max+1);
	if (entrance)
		if (strcmp(entry_type, "Undirected"))
			entrant1 = zeros(1, s_max+1);
			for (ii = 0:s_max)
				entrant1(ii+1) = sum(probE(ii+1, :).*v_guess) - 0;
			end
			FOC_E = sum(mu_guess.*entrant1);

			xE = max(GPInv(FOC_E./(w_guess*(1-subsidyE)), B_entrant), 0);
		elseif (strcmp(entry_type, "Ufuk/Sina"))
			entrant1 = zeros(1, s_max+1);
			for (ii = 0:s_max)
				entrant1(ii+1) = sum(probE(ii+1, :).*v_guess) - 0;
			end
			xE = max(GPInv(entrant1./(w_guess*(1-subsidyE)), B_entrant), 0);
		end
	end
	xE(imag(xE) ~= 0) = 0;
	
	if (entrance)
		[A, ~, b] = firm_distributions_v3_noinversion(investment_in, xE, eta, probL, probF, ...
													  prob0, ...
													  probE, prob_exog, s_max);
	else
		[~, mu_tmp, ~] = firm_distributions_v3(investment_in, xE, eta, probL, probF, prob0, ...
 										 probE, prob_exog, s_max);
		mu_guess = mu_tmp';
		
	end

	div = EIS*rho + (1-EIS)*r_guess;
	for (ii = 1:(s_max-1))
		Upsilon(s_max+1+ii) = v_guess(s_max+1+ii) - (1/div)*(profit(s_max+1+ii) - ...
				    w_guess*G(investment_in(s_max+1+ii),B)*(1-subsidy(s_max+1+ii)) + ...
				    investment_in(s_max+1+ii)*(sum(prob(s_max+1+ii, :).*v_guess) - v_guess(s_max+1+ii)) + ...
				    investment_in(s_max+1-ii)*((sum(flip(prob(s_max+1-ii, 1:s_max)).*v_guess((s_max+2):end)) + ...
										     sum(prob(s_max+1-ii,s_max+1)*v_guess(s_max+1)) + ...
										     sum(prob(s_max+1-ii, (s_max+2):end).*flip(v_guess(1:s_max)))) - ...
				    						 v_guess(s_max+1+ii)) + ...
				    eta(ii+1)*((sum(flip(prob_exog(ii, 1:s_max)).*v_guess((s_max+2):end)) + ...
								sum(prob_exog(ii, s_max+1)*v_guess(s_max+1))) - v_guess(s_max+1+ii)) + ...
				    xE(ii+1)*((sum(flip(probE(ii+1, 1:s_max)).*v_guess((s_max+2):end)) + ...
				    		   sum(probE(ii+1,s_max+1)*v_guess(s_max+1)) + ...
				    		   sum(probE(ii+1, (s_max+2):end).*flip(v_guess(1:s_max)))) - v_guess(s_max+1+ii)));

		Upsilon(s_max+1-ii) = v_guess(s_max+1-ii) - (1/div)*(profit(s_max+1-ii) - ...
				    w_guess*G(investment_in(s_max+1-ii),B)*(1-subsidy(s_max+1-ii)) + ...
				    investment_in(s_max+1-ii)*(sum(prob(s_max+1-ii, :).*v_guess) - v_guess(s_max+1-ii)) + ...
				    investment_in(s_max+1+ii)*(sum(prob(s_max+1+ii, (s_max+2):end).*flip(v_guess(1:s_max))) - ...
				    						v_guess(s_max+1-ii)) + ...
				    eta(ii+1)*(sum(prob_exog(ii,:).*v_guess) - v_guess(s_max+1-ii)) + ...
				    xE(ii+1)*(0 - v_guess(s_max+1-ii)));



	end

	Upsilon(s_max+1) = v_guess(s_max+1) - (1/div)*(profit(s_max+1) -w_guess*G(investment_in(s_max+1), B)*(1-subsidy(s_max+1)) + ...
						 investment_in(s_max+1)*(sum(prob(s_max+1, :).*v_guess) - v_guess(s_max+1)) + ...
						 investment_in(s_max+1)*(sum(prob(s_max+1, (s_max+2):end).*flip(v_guess(1:s_max))) - ...
						 					  v_guess(s_max+1)) + ...
						 xE(1)*(.5*sum(probE(1, (s_max+2):end).*flip(v_guess(1:s_max))) + .5*0 - ...
						 	    v_guess(s_max+1)));

	Upsilon(1) = v_guess(1) - (1/div)*(profit(1) -w_guess*G(investment_in(1), B)*(1-subsidy(1)) + ...
				   investment_in(1)*(sum(prob(1, :).*v_guess) - v_guess(1)) + ...
				   eta(end)*(sum(prob_exog(end,:).*v_guess) - v_guess(1)) + ...
				   xE(end)*(0 - v_guess(1)));

	Upsilon(2*s_max+1) = v_guess(end) - (1/div)*(profit(end) - w_guess*G(investment_in(end), B)*(1-subsidy(end)) + ...
					 investment_in(1)*((sum(flip(prob(1,1:s_max)).*v_guess((s_max+2):end)) + ...
					 					sum(prob(1, s_max+1)*v_guess(s_max+1)) + ...
					 					sum(prob(1, (s_max+2):end).*flip(v_guess(1,1:s_max)))) - ...
					 					v_guess(end)) + ...
					 eta(end)*((sum(flip(prob_exog(end,1:s_max)).*v_guess((s_max+2):end)) + ...
					 			sum(prob_exog(end,s_max+1)*v_guess(s_max+1))) - ...
					 			v_guess(end)) + ...
					 xE(end)*((sum(flip(probE(end,1:s_max)).*v_guess((s_max+2):end)) + ...
					 		   sum(probE(end,s_max+1)*v_guess(s_max+1)) + ...
					 		   sum(probE(end,(s_max+2):end).*flip(v_guess(1:s_max)))) - ...
					 		   v_guess(end)));


	if (entrance)
		Upsilon((2*s_max+2):(3*s_max+2)) = A*mu_guess' - b;
	end

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

	growth_rate = log(lambda)*sum(mu_guess.*(xE.*gE_iter + ...
											 flip(investment_in(1:(s_max+1))).*gF_iter + ...
											 investment_in((s_max+1):end).*gL_iter));

	if (~elastic_labor) && (EIS == 1) && (entrance)
		x0 = investment_in(s_max+1);
		xL = investment_in((s_max+2):end);
		xF = flip(investment_in(1:s_max));
		if (kappa >= 9999)
			Upsilon(3*s_max+3) = 1 - sum(mu_guess.*(G([x0, xL], B) + G([x0, xF], B) + ...
													G(xE, B_entrant) + 1./(w_guess*lambda.^min(LMS_bar_s, (0:s_max)))));
		else
			Upsilon(3*s_max+3) = 1 - sum(mu_guess.*(G([x0, xL], B) + G([x0, xF], B) + ...
													G(xE, B_entrant) + ...
													(nu_s.^(1-kappa))*(kappa-1)*(1/w_guess)./((1+nu_s.^(1-kappa)).*(kappa+nu_s.^(1-kappa))) + ...
													(nu_s.^(1-kappa))*(kappa-1)*(1/w_guess)./((1+nu_s.^(1-kappa)).*(kappa*nu_s.^(1-kappa)+1))));
			disp('CHECK THAT YOU HAVE THE RIGHT CLEARING CONDITION FOR LABOR MARKET')
		end
	elseif (elastic_labor) && (EIS ~= 1) && (entrance)
		Upsilon(3*s_max+3) = growth_rate - EIS*(r_guess - rho);
	elseif (~elastic_labor) && (EIS ~= 1) && (entrance)
		x0 = investment_in(s_max+1);
		xL = investment_in((s_max+2):end);
		xF = flip(investment_in(1:s_max));
		if (kappa >= 9999)
			Upsilon(3*s_max+3) = 1 - sum(mu_guess.*(G([x0, xL], B) + G([x0, xF], B) + ...
													G(xE, B_entrant) + 1./(w_guess*min(LMS_bar_s, lambda.^(0:s_max)))));
		else 
			Upsilon(3*s_max+3) = 1 - sum(mu_guess.*(G([x0, xL], B) + G([x0, xF], B) + ...
													G(xE, B_entrant) + ...
													(nu_s.^(1-kappa))*(kappa-1)*(1/w_guess)./((1+nu_s.^(1-kappa)).*(kappa+nu_s.^(1-kappa))) + ...
													(nu_s.^(1-kappa))*(kappa-1)*(1/w_guess)./((1+nu_s.^(1-kappa)).*(kappa*nu_s.^(1-kappa)+1))));
			disp('CHECK THAT YOU HAVE THE RIGHT CLEARING CONDITION FOR LABOR MARKET')
		end


		Upsilon(3*s_max+4) = growth_rate - EIS*(r_guess - rho);


	end




	if (~elastic_labor) && (EIS == 1) && (~entrance)
		x0 = investment_in(s_max+1);
		xL = investment_in((s_max+2):end);
		xF = flip(investment_in(1:s_max));
		if (kappa >= 9999)
			Upsilon(2*s_max+2) = 1 - sum(mu_guess.*(G([x0, xL], B) + G([x0, xF], B) + ...
													G(xE, B_entrant) + 1./(w_guess*lambda.^min(LMS_bar_s, (0:s_max)))));
		else
			Upsilon(2*s_max+2) = 1 - sum(mu_guess.*(G([x0, xL], B) + G([x0, xF], B) + ...
													G(xE, B_entrant) + ...
													(nu_s.^(1-kappa))*(kappa-1)*(1/w_guess)./((1+nu_s.^(1-kappa)).*(kappa+nu_s.^(1-kappa))) + ...
													(nu_s.^(1-kappa))*(kappa-1)*(1/w_guess)./((1+nu_s.^(1-kappa)).*(kappa*nu_s.^(1-kappa)+1))));
			disp('CHECK THAT YOU HAVE THE RIGHT CLEARING CONDITION FOR LABOR MARKET')
		end
	elseif (elastic_labor) && (EIS ~= 1) && (~entrance)
		Upsilon(2*s_max+2) = growth_rate - EIS*(r_guess - rho);
	elseif (~elastic_labor) && (EIS ~= 1) && (~entrance)
		x0 = investment_in(s_max+1);
		xL = investment_in((s_max+2):end);
		xF = flip(investment_in(1:s_max));
		if (kappa >= 9999)
			Upsilon(2*s_max+2) = 1 - sum(mu_guess.*(G([x0, xL], B) + G([x0, xF], B) + ...
													G(xE, B_entrant) + 1./(w_guess*min(LMS_bar_s, lambda.^(0:s_max)))));
		else 
			Upsilon(2*s_max+2) = 1 - sum(mu_guess.*(G([x0, xL], B) + G([x0, xF], B) + ...
													G(xE, B_entrant) + ...
													(nu_s.^(1-kappa))*(kappa-1)*(1/w_guess)./((1+nu_s.^(1-kappa)).*(kappa+nu_s.^(1-kappa))) + ...
													(nu_s.^(1-kappa))*(kappa-1)*(1/w_guess)./((1+nu_s.^(1-kappa)).*(kappa*nu_s.^(1-kappa)+1))));
			disp('CHECK THAT YOU HAVE THE RIGHT CLEARING CONDITION FOR LABOR MARKET')
		end


		Upsilon(2*s_max+3) = growth_rate - EIS*(r_guess - rho);


	end



end
