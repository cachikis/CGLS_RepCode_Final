%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: vfi_alt.m
% Author: Craig A. Chikis
% Date: 05/22/2020
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v_out = vfi_alt(value_function_in, investment_in, eta, ...
																beta, s_max, ...
																wage_share, B, gamma, xbar, entrance, ...
																subsidy, tax_rate, lambda, rho, xE, ...
																B_entrant, probL, probF, prob0, ...
																probE, prob_exog)

	

	G = @(x) (x/B).^(1/gamma);
	G_entrant = @(x) (x/B_entrant).^(1/gamma);

	xL = investment_in((s_max+2):end);
	xF = flip(investment_in(1:s_max));
	x0 = investment_in(s_max+1);

	vL = value_function_in((s_max+2):end);
	vF = flip(value_function_in(1:s_max));
	v0 = value_function_in(s_max+1);

	S = 1:s_max;

	tax_rateL = tax_rate((s_max+2):end);
	tax_rateF = flip(tax_rate(1:s_max));
	tax_rate0 = tax_rate(s_max+1);

	subsidyL = subsidy((s_max+2):end);
	subsidyF = flip(subsidy(1:s_max));
	subsidy0 = subsidy(s_max+1);

	xbar_term = 2*xbar*(~entrance) + 3*xbar*(entrance);
	eta_max = max(eta);

	profit_d = 1/(rho + xbar_term + eta_max);



	profitL = ((1 - lambda.^(-S)).*(1 - tax_rateL) - (1 - subsidyL).*G(xL)*wage_share)*profit_d;
	profitF = -(1 - subsidyF).*G(xF)*wage_share*profit_d;
	profit0 = -(1 - subsidy0)*G(x0)*wage_share*profit_d; 

	v_iter_n = 1/(xbar_term + eta_max);

	nothing0 = (1 - (2*x0 + xE(1))*v_iter_n)*v0;

	movement0 = sum(prob0*x0*v_iter_n.*value_function_in) + sum(prob0((s_max+2):end)*x0*v_iter_n.*vF) + ...
				sum(.5*probE(1, (s_max+2):end)*xE(1)*v_iter_n.*vF);
	
	nothingL = (1 - (xL + xF + eta(2:end) + xE(2:end))*v_iter_n).*vL;
	nothingF = (1 - (xL + xF + eta(2:end) + xE(2:end))*v_iter_n).*vF;

	somethingL = zeros(1, s_max);
	somethingF = zeros(1, s_max);

	no_leapfrog_master = flip(probF(:, 1:s_max), 2);
	no_leapfrog_exog_master = flip(prob_exog(:, 1:s_max),2);

	leapfrog_master = probF(:, (s_max+2):end);
	leapfrog_exog_master = prob_exog(:, (s_max+2):end);
	leader_ahead_master = probL(:, (s_max+2):end);

	tied_term_master = probF(:, s_max+1);
	tied_term_exog_master = prob_exog(:, s_max+1);

	for (s = 1:s_max)
	
		leapfrog = leapfrog_master(s,:);
		no_leapfrog = no_leapfrog_master(s,:);
		tied_term = tied_term_master(s,:);

		leapfrog_exog = leapfrog_exog_master(s,:);
		no_leapfrog_exog = no_leapfrog_exog_master(s,:);
		tied_term_exog = tied_term_exog_master(s,:);

		leader_ahead = leader_ahead_master(s,:);


		somethingL(s) = sum(probL(s,:)*xL(s)*v_iter_n.*value_function_in) + ...
									   sum((leapfrog*xF(s) + leapfrog_exog*eta(s+1))*v_iter_n.*vF + ...
										(no_leapfrog*xF(s) + no_leapfrog_exog*eta(s+1))*v_iter_n.*vL) + ...
										sum((tied_term*xF(s) + tied_term_exog*eta(s+1))*v_iter_n*v0);

		somethingF(s) = sum((probF(s,:)*xF(s) + prob_exog(s,:)*eta(s+1))*v_iter_n.*value_function_in) + ...
						sum(leader_ahead*xL(s)*v_iter_n.*vF);

		

	end


	somethingL_entrance = zeros(1, s_max);
	probE_F = flip(probE(:, 1:s_max), 2);
	probE_L = probE(:, (s_max+2):end);
	probE_0 = probE(:, s_max+1);
	for (s = 1:s_max)
		comb_noleap = probE_F(s+1, :)*xE(s+1);
		comb_leap = probE_L(s+1, :)*xE(s+1);
		comb_tied = probE_0(s+1, :)*xE(s+1);

		somethingL_entrance(s) = sum(comb_noleap*v_iter_n.*vL + comb_leap*v_iter_n.*vF) + ...
														  		 sum(comb_tied*v_iter_n*v0);
	end

	somethingL = somethingL + somethingL_entrance; 

	profitL = [profit0, profitL];
	profitF = [profit0, profitF];

	RHS_L = [movement0+nothing0, somethingL+nothingL];
	RHS_F = [movement0+nothing0, somethingF+nothingF];

	vL_out = profitL + beta.*RHS_L;
	vF_out = profitF + beta.*RHS_F;

	v_out = [flip(vF_out), vL_out(2:end)];






end




