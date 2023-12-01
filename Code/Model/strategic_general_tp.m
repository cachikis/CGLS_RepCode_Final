%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: strategic_general_tp.m
% Author: Craig A. Chikis
% Date: 06/15/2020
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = strategic_general_tp(m)


	investment = m.investment;
	investment_entrant = m.investment_entrant;
	value_functions = m.value_functions;
	eta = m.eta;
	prob = m.prob; 
	probE = m.probE;
	prob_exog = m.prob_exog; 
	s_max = m.s_max; 
	rho = m.rho; 
	wage_share = m.wage_share; 
	B = m.B;
	entrance = m.entrance; 
	subsidy = m.subsidy; 
	dx_drho_total = m.dxdrho; 
	gamma = m.gamma; 

	



	G = @(x, B) (x/B).^(1/gamma);
	GP = @(x, B) (1/(B*gamma))*(x/B).^((1/gamma)-1);
	GdP = @(x, B) (1/((B^2)*gamma))*((1/gamma)-1)*(x/B).^((1/gamma)-2);
	cs = @(x, B, tauRD) 1./(GdP(x, B)*wage_share.*(1-tauRD));


	probL = prob((s_max+2):end, :);
	probF = flip(prob(1:s_max, :), 1);
	prob0 = prob(s_max+1, :);

	probE = flip(probE, 1);
	prob_exog = flip(prob_exog, 1);


	vL = value_functions((s_max+2):end);
	vF = flip(value_functions(1:s_max));
	v0 = value_functions(s_max+1);

	xL = investment((s_max+2):end);
	xF = flip(investment(1:s_max));
	x0 = investment(s_max+1);

	xE = investment_entrant;

	eta(1) = 0;

	A = zeros((2*s_max+1)^2);
	b = zeros((2*s_max+1)^2, 1);


	follower_idx_v = 1:(s_max*(2*s_max+1));
	Sigma_follower = repelem(-s_max:s_max, s_max);
	S_follower = repmat(flip(1:s_max), 1, 2*s_max+1);

	tied_idx_v = repelem((1:(2*s_max+1)) + follower_idx_v(end), 1, s_max);

	leader_idx_v = (1:(s_max*(2*s_max+1))) + tied_idx_v(end);
	Sigma_leader = repelem(-s_max:s_max, s_max);
	S_leader = repmat(1:s_max, 1, 2*s_max + 1);

	Sigma_tied = -s_max:s_max;


	for (i = 1:length(follower_idx_v))

		A(follower_idx_v(i), follower_idx_v(i)) = rho + xL(S_follower(i)) + xF(S_follower(i)) + ...
												  xE(S_follower(i)+1) + eta(S_follower(i)+1);


		L_term = probL(S_follower(i), :)*xL(S_follower(i));
		F_term = probF(S_follower(i), :)*xF(S_follower(i));
		exog_term = prob_exog(S_follower(i), :)*eta(S_follower(i)+1);

		L_cost = flip(L_term((s_max+2):end));
		F_F = F_term(1:s_max);
		F_0 = F_term(s_max+1);
		F_lp = F_term((s_max+2):end);
		exog_F = exog_term(1:s_max);
		exog_0 = exog_term(s_max+1);
		


		A(follower_idx_v(i), find(Sigma_follower == Sigma_follower(i))) =  ...
				A(follower_idx_v(i), find(Sigma_follower == Sigma_follower(i))) - (L_cost + F_F + exog_F);


		A(follower_idx_v(i), (2*s_max+1)*(s_max+1)+find(Sigma_follower == Sigma_follower(i))) = ...
				A(follower_idx_v(i), (2*s_max+1)*(s_max+1)+find(Sigma_follower == Sigma_follower(i))) - F_lp;

		A(follower_idx_v(i), tied_idx_v(i)) =  ... 
				A(follower_idx_v(i), tied_idx_v(i)) - (exog_0 + F_0);


		b_term = sum(flip(probL(S_follower(i), 1:s_max)).*vL) + sum(probL(S_follower(i), (s_max+2):end).*vF) + ...
				 	probL(S_follower(i), s_max+1)*v0 - vF(S_follower(i));



		b(follower_idx_v(i)) = (Sigma_follower(i) == S_follower(i))*b_term;

	end

	for (i = 1:length(leader_idx_v))

		A(leader_idx_v(i), leader_idx_v(i)) = rho + xL(S_leader(i)) + xF(S_leader(i)) + ...
												  xE(S_leader(i)+1) + eta(S_leader(i)+1);


		L_term = probL(S_leader(i), :)*xL(S_leader(i));
		F_term = probF(S_leader(i), :)*xF(S_leader(i));
		exog_term = prob_exog(S_leader(i), :)*eta(S_leader(i)+1);
		E_term = probE(S_leader(i)+1, :)*xE(S_leader(i)+1);

		L_gain = L_term((s_max+2):end);
		F_F = flip(F_term(1:s_max));
		F_0 = F_term(s_max+1);
		F_lp = F_term((s_max+2):end);
		exog_F = flip(exog_term(1:s_max));
		exog_0 = exog_term(s_max+1);
		E_F = flip(E_term(1:s_max));
		E_lp = E_term((s_max+2):end);
		E_0 = E_term(s_max+1);


		A(leader_idx_v(i), (2*s_max+1)*(s_max+1) + find(Sigma_leader == Sigma_leader(i))) =  ...
				A(leader_idx_v(i), (2*s_max+1)*(s_max+1)+find(Sigma_follower == Sigma_follower(i))) - (L_gain + F_F + exog_F + E_F);


		A(leader_idx_v(i), find(Sigma_leader == Sigma_leader(i))) = ...
				A(leader_idx_v(i), find(Sigma_leader == Sigma_leader(i))) - (F_lp + E_lp);

		A(leader_idx_v(i), tied_idx_v(i)) =  ... 
				A(leader_idx_v(i), tied_idx_v(i)) - (exog_0 + F_0 + E_0);


		b_term = sum(flip(probF(S_leader(i), 1:s_max)).*vL) + sum(probF(S_leader(i), (s_max+2):end).*vF) + ...
				 (probL(S_leader(i), s_max+1))*v0 - vL(S_leader(i));



		b(leader_idx_v(i)) = (Sigma_leader(i) == -S_leader(i))*b_term;

	end

	tied_idx_v_unique = unique(tied_idx_v);

	for (i = 1:length(tied_idx_v_unique))

		A(tied_idx_v_unique(i), tied_idx_v_unique(i)) = rho + 2*x0 + xE(1); 
												 
		T_term = prob0*x0;
		T_c_term = prob0*x0;
		E_term = probE(1, :)*xE(1);

		T_gain = T_term((s_max+2):end);
		T_c_loss = flip(T_c_term((s_max+2):end));
		E_loss = flip(.5*E_term((s_max+2):end));

		A(tied_idx_v_unique(i), (2*s_max+1)*(s_max+1)+find(Sigma_leader == Sigma_tied(i))) = ...
				A(tied_idx_v_unique(i), (2*s_max+1)*(s_max+1)+find(Sigma_leader == Sigma_leader(i))) - T_gain;

		A(tied_idx_v_unique(i), find(Sigma_follower == Sigma_tied(i))) = ...
				A(tied_idx_v_unique(i), find(Sigma_follower == Sigma_follower(i))) - (T_c_loss + E_loss);

		b_term = sum(flip(prob0(1:s_max)).*vF) + sum(prob0((s_max+2):end).*vF)  - v0;

		b(tied_idx_v_unique(i)) = (Sigma_tied(i) == 0)*b_term;	


	
	end
	A_sparse = sparse(A);
	parv_parxsig = A_sparse \ b;

	parv_F = reshape(parv_parxsig(1:follower_idx_v(end)), s_max, 2*s_max+1);
	parv_L = reshape(parv_parxsig(leader_idx_v(1):leader_idx_v(end)), s_max, 2*s_max+1);
	parv_0 = reshape(parv_parxsig(tied_idx_v_unique(1):tied_idx_v_unique(end)), 2*s_max+1, 1);

	parv_F = parv_F';
	parv_L = parv_L';

	parv_all = [parv_F, parv_0, parv_L]; 


	parDv_F = zeros(size(parv_F));
	parDv_L = zeros(size(parv_L));
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
		strategic_summation(s) = dx_drho_total*parx_parxsig(:, s); 
	end

	m.strategic_summation = strategic_summation; 
	m.parv_parxsig = parv_all; 
	m.parDv_parxsig = parDv_all;
	m.parx_parxsig = parx_parxsig; 




end
