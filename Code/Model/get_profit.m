%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: get_profit.m
% Author: Craig A. Chikis
% Note(s):
% MA-MFA
% Date: 11/24/20202
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [profit, nu_s] = get_profit(s_max, lambda, kappa, tax_rate, LMS_bar_s, two_adjust)


	options_fzero = optimset('Display', 'off');


	nu_s = nan(1, s_max+1);
	fval = nan(1, s_max+1);
	flag = nan(1, s_max+1);
	for (s = 0:s_max) 
		eq_solve = @(rho) solve_rho_impl(rho, lambda, s, kappa);
		[nu_s(s+1), fval(s+1), flag(s+1)] = fzero(eq_solve, 1, options_fzero);
	end

	if (ismember(-1, sign(flag)))
		disp('Probem! Negative exit flag from fsolve from nu_s solution!')
		return;
	end

	if (kappa >= 9999) 
		profitL = (1 - lambda.^(-(1:s_max))).*(1 - tax_rate((s_max+2):end));
		profitF = zeros(1,s_max);
		profit0 = 0;

		profitL(:) = profitL(min((1:s_max)+0, LMS_bar_s));
		profitF(:) = profitF(min((1:s_max)+0, LMS_bar_s));
	else
		profitL = (1-tax_rate((s_max+1):end)).*(nu_s.^(1-kappa))./(kappa + nu_s.^(1-kappa));
		profitF = (1-flip(tax_rate(1:(s_max+1))))./(kappa*nu_s.^(1-kappa) + 1);
		profit0 = (1-tax_rate(s_max+1))/(kappa+1);

		profitL(:) = profitL(min((0:s_max)+1, LMS_bar_s+1));
		profitF(:) = profitF(min((0:s_max)+1, LMS_bar_s+1));

		profitL = profitL(2:end);
		profitF = profitF(2:end);
	end

	profit = two_adjust*[flip(profitF), profit0, profitL];
	nu_s(:) = nu_s(min((0:s_max)+1, LMS_bar_s+1));



end
