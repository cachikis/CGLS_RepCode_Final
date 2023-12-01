%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: entrant_strategic.m
% Author: Craig A. Chikis
% Date: 06/23/2020
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = incumbent_entrant_strategic(m)

investment = m.investment; 
investment_entrant = m.investment_entrant; 
eta = m.eta;
rho = m.rho; 
prob = m.prob;  
probE = m.probE; 
prob_exog = m.prob_exog; 
entry_type = m.entry_type; 
s_max = m.s_max; 
value_functions = m.value_functions; 
wage_share = m.wage_share; 
subsidy = m.subsidy; 
B = m.B;
dx_drho_total = m.dxdrho;  
gamma = m.gamma; 



if (strcmp(entry_type, "Undirected"))
	G = @(x, B) (x/B).^(1/gamma);
	GP = @(x, B) (1/(B*gamma))*(x/B).^((1/gamma)-1);
	GdP = @(x, B) (1/((B^2)*gamma))*((1/gamma)-1)*(x/B).^((1/gamma)-2);
	cs = @(x, B, tauRD) 1./(GdP(x, B)*wage_share.*(1-tauRD));

	A = zeros(2*s_max+1);
	b = zeros(size(A,1), 1);

	xL = investment((s_max+2):end);
	xF = flip(investment(1:s_max));
	x0 = investment(s_max+1);

	vL = value_functions((s_max+2):end);
	vF = flip(value_functions(1:s_max));
	v0 = value_functions(s_max+1);

	probL = prob((s_max+2):end, :);
	probF = flip(prob(1:s_max, :), 1);
	prob0 = prob(s_max+1, :);

	probE = flip(probE, 1);
	prob_exog = flip(prob_exog, 1);

	xE = investment_entrant;

	for (ii = 1:length(1:s_max))
		A(ii, ii) = rho + xL(s_max+1-ii) + xF(s_max+1-ii) + eta(s_max+2-ii) + xE(s_max+2-ii);

		L_term = probL(s_max+1-ii, :)*xL(s_max+1-ii);
		F_term = probF(s_max+1-ii, :)*xF(s_max+1-ii);
		exog_term = prob_exog(s_max+1-ii, :)*eta(s_max+2-ii);

		L_cost = flip(L_term((s_max+2):end));
		F_F = F_term(1:s_max);
		F_0 = F_term(s_max+1);
		F_lp = F_term((s_max+2):end);
		exog_F = exog_term(1:s_max);
		exog_0 = exog_term(s_max+1);


		A(ii, 1:s_max) = A(ii, 1:s_max) + -(L_cost + F_F + exog_F);
		A(ii, (s_max+2):end) = A(ii, (s_max+2):end) + -F_lp;
		A(ii, s_max+1) = A(ii, s_max+1) + -(exog_0 + F_0);

		b_term = -value_functions(ii);

		b(ii) = b_term;

	end

	for (ii = 1:length(1:s_max))
		A(s_max+1+ii, s_max+1+ii) = rho + xL(ii) + xF(ii) + eta(ii+1) + xE(ii+1);

		L_term = probL(ii, :)*xL(ii);
		F_term = probF(ii, :)*xF(ii);
		exog_term = prob_exog(ii, :)*eta(ii+1);
		E_term = probE(ii+1, :)*xE(ii+1);

		L_gain = L_term((s_max+2):end);
		F_F = flip(F_term(1:s_max));
		F_0 = F_term(s_max+1);
		F_lp = F_term((s_max+2):end);
		exog_F = flip(exog_term(1:s_max));
		exog_0 = exog_term(s_max+1);
		E_F = flip(E_term(1:s_max));
		E_lp = E_term((s_max+2):end);
		E_0 = E_term(s_max+1);

		A(ii+s_max+1, (s_max+2):end) = A(ii+s_max+1, (s_max+2):end) + -(L_gain + F_F + exog_F + E_F);
		A(ii+s_max+1, 1:s_max) = A(ii+s_max+1, 1:s_max) + -(F_lp + E_lp);
		A(ii+s_max+1, s_max+1) = A(ii+s_max+1, s_max+1) + -(exog_0 + F_0 + E_0);

		b_term = sum(flip(probE(ii+1, 1:s_max)).*vL) + sum(probE(ii+1, (s_max+2):end).*vF) + ...
				 (probE(ii+1, s_max+1))*v0 - vL(ii);

		b(s_max+1+ii) = b_term;


	end


	A(s_max+1, s_max+1) = rho + 2*x0 + xE(1);
	T_term = prob0*x0;
	T_c_term = prob0*x0;
	E_term = probE(1,:)*xE(1);


	T_gain = T_term((s_max+2):end);
	T_c_loss = flip(T_c_term((s_max+2):end));
	E_loss = flip(.5*E_term((s_max+2):end));

	A(s_max+1, (s_max+2):end) = -T_gain;
	A(s_max+1, 1:s_max) = -(T_c_loss + E_loss);

	b_term = .5*(sum(flip(probE(1, 1:s_max)).*vF) + sum(probE(1, (s_max+2):end).*vF) + probE(1, s_max+1)*v0)  - v0;

	b(s_max+1) = b_term;

	parvs_parxe = A \ b;

	parv_F = parvs_parxe(1:s_max);
	parv_0 = parvs_parxe(s_max+1);
	parv_L = parvs_parxe((s_max+2):end);

	parv_all = parvs_parxe';

	parDv_F = zeros(size(parv_F))';
	parDv_L = zeros(size(parv_L))';
	for (s = 1:s_max)
		frac_itF = prob(s_max+1-s, :);
		frac_itL = prob(s_max+1+s, :);

		pc_f = 0;
		pc_l = 0;
		for (f = 1:length(frac_itF))
			pc_f = pc_f + frac_itF(f)*parv_all(:, f);
			pc_l = pc_l + frac_itL(f)*parv_all(:, f);
		end

		parDv_F(:, s) = pc_f - parv_all(:, s_max+1-s);
		parDv_L(:, s) = pc_l - parv_all(:, s_max+1+s);
	end

	frac_it0 = prob(s_max+1, :);
	pc_0 = 0;
	for (f = 1:length(frac_it0))
		pc_0 = pc_0 + frac_it0(f)*parv_all(:, f);
	end

	parDv_0 = pc_0 - parv_all(:, s_max+1);

	parDv_all = [flip(parDv_F, 2), parDv_0, parDv_L];

	parx_parxsig = zeros(size(parDv_all));
	for (s = 1:size(parDv_all,2))
		parx_parxsig(:, s) = cs(investment(s), B, subsidy(s))*parDv_all(:, s);
	end

	for (s = 1:size(parx_parxsig, 2))
		strategic_summation(s) = dx_drho_total(1)*parx_parxsig(:, s);
	end

else

	A = zeros((s_max+1)*(2*s_max+1));
	b = zeros(size(A,1), 1);

	parv_all = zeros(1, 2*s_max+1);
	A = zeros(2*s_max+1);
	b = zeros(2*s_max+1,1); 
	strategic_summation = zeros(1,2*s_max+1); 
	parx_parxsig = zeros(1,2*s_max+1);




end


m.parv_parxe = parv_all; 
m.parx_parxesig = parx_parxsig;



end
