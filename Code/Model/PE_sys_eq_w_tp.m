%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: PE_sys_eq_w_tp.m
% Author: Craig A. Chikis
% Date: 06/11/2020
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = PE_sys_eq_w_tp(m)


	value_functions = m.value_functions;
	investment = m.investment;
	investment_entrant = m.investment_entrant; 
	s_max = m.s_max; 
	rho = m.rho; 
	eta = m.eta;
	wage_share = m.wage_share; 
	prob = m.prob; 
	probE = m.probE; 
	prob_exog = m.prob_exog; 
	B = m.B;
	gamma = m.gamma; 
	xbar = m.xbar; 
	lambda = m.lambda;
	subsidy = m.subsidy; 

	G = @(x, B) (x/B).^(1/gamma);
	GP = @(x, B) (1/(B*gamma))*(x/B).^((1/gamma)-1);
	GdP = @(x, B) (1/((B^2)*gamma))*((1/gamma)-1)*(x/B).^((1/gamma)-2);
	cs = @(x, B, subsidy) 1./(GdP(x, B)*wage_share.*(1-subsidy));


	vL = value_functions((s_max+2):end);
	vF = flip(value_functions(1:s_max));
	v0 = value_functions(s_max+1);

	xL = investment((s_max+2):end);
	xF = flip(investment(1:s_max));
	x0 = investment(s_max+1);

	xE = investment_entrant;

	A_parv_parrho = zeros(2*s_max+1);
	b_parv_parrho = zeros(2*s_max+1, 1);

	b_parv_parw = zeros(2*s_max+1, 1);

	probL = prob((s_max+2):end, :);
	probF = flip(prob(1:s_max, :), 1);
	prob0 = prob(s_max+1, :);

	prob_exog = flip(prob_exog, 1);
	probE = flip(probE, 1);
	

	for (s = 1:s_max)

		probF_it = probF(s, :);
		probL_it = probL(s, :);
		probE_it = probE(s+1, :); 
		prob_exog_it = prob_exog(s, :);

		F_iter_forF = probF_it*xF(s) + prob_exog_it*eta(s+1);
		L_iter_forF = flip(probL_it((s_max+2):end)*xL(s)); 

		A_parv_parrho(s_max+1-s, s_max+1-s) = A_parv_parrho(s_max+1-s, s_max+1-s) + ...
												rho + xL(s) + xF(s) + xE(s+1) + eta(s+1);

		
		A_parv_parrho(s_max+1-s, 1:s_max) = A_parv_parrho(s_max+1-s, 1:s_max) - L_iter_forF;
		A_parv_parrho(s_max+1-s, :) = A_parv_parrho(s_max+1-s, :) - F_iter_forF;

		b_parv_parrho(s_max+1-s) = -vF(s);
		b_parv_parw(s_max+1-s) = -(1-subsidy(s_max+1-s))*G(xF(s), B);


		A_parv_parrho(s_max+1+s, s_max+1+s) = A_parv_parrho(s_max+1+s, s_max+1+s) + ...
												rho + xL(s) + xF(s) + xE(s+1) + eta(s+1);

		F_iter_forL = flip(probF_it*xF(s) + prob_exog_it*eta(s+1) + probE_it*xE(s+1));
		L_iter_forL = probL_it((s_max+2):end)*xL(s);

		A_parv_parrho(s_max+1+s, :) = A_parv_parrho(s_max+1+s, :) - F_iter_forL;
		A_parv_parrho(s_max+1+s, (s_max+2):end) = A_parv_parrho(s_max+1+s, (s_max+2):end) - L_iter_forL;

		b_parv_parrho(s_max+1+s) = -vL(s);
		b_parv_parw(s_max+1+s) = -(1-subsidy(s_max+1+s))*G(xL(s), B);




	end

	b_parv_parrho(s_max+1) = -v0;
	b_parv_parw(s_max+1) = -(1-subsidy(s_max+1))*G(x0, B);

	A_parv_parrho(s_max+1, s_max+1) = rho + 2*x0 + xE(1);

	F_iter_for0 = flip(prob0((s_max+2):end))*x0 + flip(probE(1,(s_max+2):end))*.5*xE(1);
	L_iter_forL = prob0((s_max+2):end)*x0;

	A_parv_parrho(s_max+1, 1:s_max) = A_parv_parrho(s_max+1, 1:s_max) - F_iter_for0;
	A_parv_parrho(s_max+1, (s_max+2):end) = A_parv_parrho(s_max+1, (s_max+2):end) - L_iter_forL;


	A_parv_parw = A_parv_parrho;
	par_Dv_parrho = (A_parv_parrho \ b_parv_parrho)';
	par_Dv_parw = (A_parv_parw \ b_parv_parw)';

	FOC_F = zeros(1, length(1:s_max));
	FOC_L = zeros(1, length(1:s_max));
	FOC_L2 = zeros(size(FOC_L));
	FOC_F2 = zeros(size(FOC_F));
	for (s = 1:s_max)
		FOC_F(s) = sum(probF(s, :).*par_Dv_parrho) - par_Dv_parrho(s_max+1-s);
		FOC_L(s) = sum(probL(s, :).*par_Dv_parrho) - par_Dv_parrho(s_max+1+s);

		FOC_L2(s) = sum(probL(s, :).*par_Dv_parw) - par_Dv_parw(s_max+1+s);
		FOC_F2(s) = sum(probF(s, :).*par_Dv_parw)- par_Dv_parw(s_max+1-s);

	end
	FOC_0 = sum(prob0.*par_Dv_parrho) - par_Dv_parrho(s_max+1);
	FOC_02 = sum(prob0.*par_Dv_parw) - par_Dv_parw(s_max+1);

	parx_parrho = cs(investment, B, subsidy).*[flip(FOC_F), FOC_0, FOC_L];
	parx_parw = cs(investment, B, subsidy).*(-GP(investment, B).*(1-subsidy) + [flip(FOC_F2), FOC_02, FOC_L2]);


	parrho_FOC = [flip(FOC_F), FOC_0, FOC_L];
	parw_FOC = [flip(FOC_F2), FOC_02, FOC_L2];


	m.parx_parrho = parx_parrho;
	m.parx_parw = parx_parw;
	m.par_Dv_parrho = par_Dv_parrho; 
	m.par_Dv_parw = par_Dv_parw; 
	m.parrho_FOC = parrho_FOC;
	m.parw_FOC = parw_FOC; 


end
