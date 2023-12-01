%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: get_BGP.m
% Author: Craig A. Chikis
% Date: 11/30/2020
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [value_functions, investment, firm_distributions, ...
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
		  										entrance, entry_type, subsidyE, kappa, alpha_fric,EIS, ...
		  										LMS_bar_s, eta)

	
	s_max = (size(prob,2)-1)/2;

	probE_iter = flip(probE, 1);
	prob_exog_iter = flip(prob_exog, 1);
	prob0 = prob(s_max+1,:);

	G = @(x,B) (x/B).^(1/gamma);
	GPInv = @(z,B) B*((z*B*gamma).^(gamma/(1-gamma)));
	GInv = @(z,B) B*z.^gamma;

	value_functions = v_in(1:(2*s_max+1));
	if (entrance)
		firm_distributions = v_in((2*s_max+2):(3*s_max+2));
	end

	
	if (elastic_labor) && (EIS == 1) && (entrance)
		wage_share = 1;
		real_interest_rate = nan;
	elseif (elastic_labor) && (EIS ~= 1) && (entrance)
		wage_share = 1;
		real_interest_rate = v_in(3*s_max+3);
	elseif (~elastic_labor) && (EIS == 1) && (entrance)
		wage_share = v_in(3*s_max+3);
		real_interest_rate = nan;
	elseif (~elastic_labor) && (EIS ~= 1) && (entrance)
		wage_share = v_in(3*s_max+3);
		real_interest_rate = v_in(3*s_max+4);
	end

	if (elastic_labor) && (EIS == 1) && (~entrance)
		wage_share = 1;
		real_interest_rate = nan;
	elseif (elastic_labor) && (EIS ~= 1) && (~entrance)
		wage_share = 1;
		real_interest_rate = v_in(2*s_max+2);
	elseif (~elastic_labor) && (EIS == 1) && (~entrance)
		wage_share = v_in(2*s_max+2);
		real_interest_rate = nan;
	elseif (~elastic_labor) && (EIS ~= 1) && (~entrance)
		wage_share = v_in(2*s_max+2);
		real_interest_rate = v_in(2*s_max+3);
	end


	investment = zeros(1, 2*s_max+1);
	FOC = zeros(1, 2*s_max+1);
	for (ii = 1:s_max)
		FOC(s_max+1+ii) = sum(value_functions.*prob(s_max+1+ii, :)) - value_functions(s_max+1+ii);
		FOC(s_max+1-ii) = sum(value_functions.*prob(s_max+1-ii, :)) - value_functions(s_max+1-ii);
	end
	FOC(s_max+1) = sum(value_functions.*prob(s_max+1,:)) - value_functions(s_max+1);
	investment = max(GPInv(FOC./(wage_share*(1-subsidy)), B), 0);
	investment(imag(investment) ~= 0) = 0;


	finout = GInv(alpha_fric*value_functions./(wage_share*(1-subsidy)), B);
	min_out = min(investment, finout); 

	fric_binding = (min_out < investment); 

	investment = max(min_out, 0);

	entrant1 = zeros(1,s_max+1);
	for (ii = 0:s_max)
		entrant1(ii+1) = sum(value_functions.*probE_iter(ii+1,:)) - 0;
	end

	xE = zeros(1, s_max+1);
	if (entrance)
		if (strcmp(entry_type, "Ufuk/Sina"))
			FOC_E = entrant1;
			xE = max(GPInv(FOC_E./(wage_share*(1-subsidyE)), B_entrant), 0);
		elseif (strcmp(entry_type, "Undirected"))
			FOC_E = sum(firm_distributions.*entrant1);
			xE = max(GPInv(FOC_E/(wage_share*(1-subsidy(2))), B_entrant), 0);
			xE = repelem(xE, length(0:s_max));
		end
	end
	xE(imag(xE) ~= 0) = 0;
	investment_entrant = xE;
	
	if (~entrance)
		probL = prob((s_max+2):end, :);
		probF = flip(prob(1:s_max, :), 1);
		prob0 = prob(s_max+1, :);
		[~, mu_tmp, ~] = firm_distributions_v3(investment, investment_entrant, eta, probL, probF, prob0, ...
 										   probE_iter, prob_exog_iter, s_max);
		firm_distributions = mu_tmp';
	end



	gE_iter = zeros(1, length(0:s_max));
	gL_iter = zeros(1, length(0:s_max));
	gF_iter = zeros(1, length(0:s_max));

	probL_iter = prob((s_max+1):end,:);
	probF_iter = flip(prob(1:(s_max+1), :),1);
		

	for (s = 0:s_max)
		gE_iter(s+1) = sum(probE_iter(s+1,(s_max+2):end).*(1:s_max));
		gL_iter(s+1) = sum(probL_iter(s+1,(s_max+2+s):end).*(((1+s):s_max)-s));
		gF_iter(s+1) = sum(probF_iter(s+1, (s_max+2):end).*(1:s_max));
	end

	growth_rate_LMS = log(lambda)*sum(firm_distributions.*(investment_entrant.*gE_iter.*((0:s_max) <= LMS_bar_s) + ...
													   flip(investment(1:(s_max+1))).*gF_iter.*((0:s_max) <= LMS_bar_s) + ...
													   investment((s_max+1):end).*gL_iter.*((0:s_max) <= LMS_bar_s)));

	growth_rate = log(lambda)*sum(firm_distributions.*(investment_entrant.*gE_iter + ...
													   flip(investment(1:(s_max+1))).*gF_iter + ...
													   investment((s_max+1):end).*gL_iter));

	markup_array = lambda.^(0:s_max);
	markup = sum(firm_distributions.*markup_array) - 1;

	markup_eq = nan;
	markup_eq_array = nan(size(markup_array));

	markup_wt = nan;
	markup_wt_array = nan(1,s_max+1);

	lab_demand = cell(1,3);


	if (~(kappa >= 9999)) && (elastic_labor)
		del_s = nu_s.^(1-kappa)./(nu_s.^(1-kappa) + 1);
		del_ms = 1./(nu_s.^(1-kappa)+1);

		m_s = (kappa + nu_s.^(1-kappa))/(kappa-1);
		m_ms = (kappa*nu_s.^(1-kappa) + 1)./((kappa-1)*nu_s.^(1-kappa));

		markup_eq_array = .5*m_s + .5*m_ms;
		markup_wt_array = del_s.*m_s + del_ms.*m_ms;

		markup_eq = sum(firm_distributions.*markup_eq_array);
		markup_wt = sum(firm_distributions.*markup_wt_array);

		lab_demand{1} = (1-subsidy).*[G(investment(1:s_max), B), G(investment((s_max+1):end),B)];
		lab_demand{2} = [flip(del_ms(2:end) - 1./(kappa*nu_s(2:end).^(1-kappa) + 1)), ...
				     	 del_s - nu_s.^(1-kappa)./(kappa + nu_s.^(1-kappa))];
		lab_demand{3} = (1-subsidyE).*G(investment_entrant, B_entrant);
	elseif (kappa >= 9999) 
		lab_demand{1} = (1-subsidy).*[G(investment(1:s_max),B), G(investment((s_max+1):end), B)];
		lab_demand{2} = [zeros(1,s_max), .5/(wage_share*lambda^min(0, LMS_bar_s)), ...
							1./(wage_share*lambda.^min(LMS_bar_s, (1:s_max)))];
		lab_demand{3} = (1-subsidyE).*G(investment_entrant, B_entrant);
	end

	

	RMSE = nan;
	value_flag = nan;
	labor_flag = nan;

	

	if (isnan(real_interest_rate))
		real_interest_rate = rho + growth_rate/EIS;
	end
	




end
