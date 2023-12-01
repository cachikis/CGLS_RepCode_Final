%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: AA2012_scu.m
% Author: Craig A. Chikis
% Date: 01/11/2021
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = AA2012_scu(s_max, phi_wt_exog, eta)


	two_adjust = 1;
	w_guess = 1;
	v_guess = ones(1,2*s_max+1);
	rho = power(1.05,1/12)-1;
	xbar = 1;
	lab_tol = 1e-4;
	val_tol = 5e-8;
	B = 0.1;
	gamma = 0.35;
	lambda = 1.05;
	entrance = false;
	subsidy = repelem(0, 2*s_max+1);
	tax_rate = repelem(0, 2*s_max+1);
	subsidyE = subsidy(1:(s_max+1));
	B_entrant = 0.050514;
	PE = false;
	entry_type = "Undirected";
	model_str = "entry_current";
	phi_wt = 0;
	phi_wt_tilde = 0.479114;
	l = 0;
	l_tilde = 0;
	l_opt = 0;
	l_exog = 0;
	kappa = 9999;
	LMS_bar_s = Inf;
	elastic_labor = false;
	algorithm_fsolve = 'trust-region-dogleg';
	alpha_fric = Inf;
	bridge_method = false;
	EIS = 1;
	r_guess = power(1.03,1/12)-1;
	winsor_vec = [1, 99];

	[prob, probE, prob_exog] = F_mat_main(phi_wt, phi_wt_tilde, phi_wt_exog, l, l_tilde, s_max);

	m.s_max = s_max;
	m.w_guess = w_guess;
	m.v_guess = v_guess;
	m.rho = rho;
	m.xbar = xbar;
	m.eta = eta;
	m.lab_tol = lab_tol;
	m.val_tol = val_tol; 
	m.B = B;
	m.gamma = gamma; 
	m.lambda = lambda;  
	m.entrance = entrance; 
	m.subsidy = subsidy; 
	m.tax_rate = tax_rate;
	m.subsidyE = subsidyE;
	m.B_entrant = B_entrant;
	m.PE = PE;
	m.entry_type = entry_type;
	m.prob = prob; 
	m.probE = probE; 
	m.prob_exog = prob_exog; 
	m.phi_wt = phi_wt;
	m.phi_wt_exog = phi_wt_exog;
	m.phi_wt_tilde = phi_wt_tilde; 
	m.l = l;
	m.l_tilde = 0;
	m.kappa = kappa; 
	m.LMS_bar_s = LMS_bar_s; 
	m.elastic_labor = elastic_labor; 
	m.algorithm_fsolve = algorithm_fsolve; 
	m.alpha_fric = alpha_fric; 
	m.bridge_method = bridge_method; 
	m.EIS = EIS; 
 	m.r_guess = r_guess; 
	m.v_max_try = 100000; 
	m.two_adjust = two_adjust; 







end
