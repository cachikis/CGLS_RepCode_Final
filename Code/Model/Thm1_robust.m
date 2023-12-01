%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: Thm1_robust.m
% Author: Craig A. Chikis
% Date: 12/08/2020
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = Thm1_robust(m, block)



firm_distributions = m.firm_distributions;
investment = m.investment;
investment_entrant = m.investment_entrant;
wage_share = m.wage_share; 
lambda = m.lambda; 
s_max = m.s_max;
B = m.B;
gamma = m.gamma;
B_entrant = m.B_entrant;
subsidy = m.subsidy;
prob = m.prob;
probE = m.probE;
prob_exog = m.prob_exog; 
parv_paromega = m.par_Dv_parw;
Dv_omega = m.Dv_omega; 
subsidyE = m.subsidyE; 
parv_parxe = m.parv_parxe; 
parv_parxsig = m.parv_parxsig; 
eta = m.eta;
value_functions = m.value_functions;
entrance = m.entrance; 
entry_type = m.entry_type;
elastic_labor = m.elastic_labor; 
LMS_bar_s = m.LMS_bar_s; 

eta(1) = 0;
	
	probL = prob((s_max+2):end, :);
	probF = flip(prob(1:s_max, :), 1);
	prob0 = prob(s_max+1, :);

	probE = flip(probE, 1);
	prob_exog = flip(prob_exog, 1);

	G = @(x, B) (x/B).^(1/gamma);
	GP = @(x, B) (1/(B*gamma))*(x/B).^((1/gamma)-1);
	GdP = @(x, B) (1/((B^2)*gamma))*((1/gamma)-1)*(x/B).^((1/gamma)-2);
	cs = @(x, B, subsidy) 1./(GdP(x, B)*wage_share.*(1-subsidy));

	if (strcmp(entry_type, "Ufuk/Sina"))
		H_Upsilon = zeros(4*s_max+5);
	else
		H_Upsilon = zeros(3*s_max+5);
	end


	H_Upsilon(1, 2) = (-1/(wage_share^2))*sum(firm_distributions.*(lambda.^(-min(LMS_bar_s, (0:s_max)))));
	H_Upsilon(1, 3:(s_max+2)) = flip(firm_distributions(2:end)).*GP(investment(1:s_max), B);
	H_Upsilon(1, s_max+3) = 2*firm_distributions(1)*GP(investment(s_max+1), B);
	H_Upsilon(1, (s_max+4):(2*s_max+3)) = firm_distributions(2:end).*GP(investment((s_max+2):end), B);
	H_Upsilon(1, (2*s_max+4)) = GP(investment_entrant(1), B_entrant);
	H_Upsilon(1, 2*s_max+5) = (1/(wage_share)) + 2*G(investment(s_max+1),B);
	H_Upsilon(1, (2*s_max+6):end) = (1./(wage_share*(lambda.^min(LMS_bar_s, (1:s_max))))) + G(investment((s_max+2):end), B) + ...
									G(flip(investment(1:s_max)), B);


	H_Upsilon(2,1) = 1;

	for (ii = 1:s_max)
		H_Upsilon(2, 2+ii) = -log(lambda)*firm_distributions(end+1-ii)*sum(probF(end+1-ii, (s_max+2):end).*(1:s_max));
		H_Upsilon(2, s_max+3+ii) = -log(lambda)*firm_distributions(1+ii)*sum(probL(ii, (s_max+2+ii):end).*(((1+ii):s_max)-ii));
	end

	H_Upsilon(2, s_max+3) = -2*log(lambda)*firm_distributions(1)*sum(prob0((s_max+2):end).*(1:s_max));


	other_term = zeros(1, s_max+1);

	for (ii = 1:length(other_term))
		other_term(ii) = sum(probE(ii, (s_max+2):end).*(1:s_max)); 
	end

	H_Upsilon(2, 2*s_max+4) = -log(lambda)*sum(firm_distributions.*other_term);



	for (ii = 0:s_max)
		H_Upsilon(2, 2*s_max+5+ii) = -log(lambda)*(sum((investment_entrant(ii+1)*probE(ii+1, (s_max+2):end) + ...
														investment(s_max+1-ii)*prob(s_max+1-ii, (s_max+2):end)).*(1:s_max)) + ...
												   investment(s_max+1+ii)*sum(prob(s_max+1+ii, (s_max+2+ii):end).*(((ii+1):s_max) - ii)));
	end
	H_Upsilon(3, 3:(s_max+2)) = prob(1:s_max, s_max+1)'.*flip(firm_distributions(2:end));
	H_Upsilon(3, s_max+3) = -2*firm_distributions(1);
	H_Upsilon(3, (s_max+4):(2*s_max+3)) = 0;
	H_Upsilon(3, 2*s_max+4) = ...
		-firm_distributions(1) + sum(probE(2:end, s_max+1)'.*firm_distributions(2:end)); 

	H_Upsilon(3, 2*s_max+5) = -(2*investment(s_max+1) + investment_entrant(1));
	H_Upsilon(3, (2*s_max+6):end) = investment((s_max+2):end).*prob((s_max+2):end, s_max+1)' + ...
									investment_entrant(2:end).*probE(2:end, s_max+1)' + ...
									flip(investment(1:s_max)).*flip(prob(1:s_max, s_max+1),1)' + ...
									eta(2:end).*prob_exog(:, s_max+1)';

	

	for (s = 1:(s_max-1))
		rownum = s+3;

		enter_Ls = [prob((s_max+1):end, s_max+1+s), prob((s_max+1):end, s_max+1-s)];
		enter_Fs = [flip(prob(1:(s_max+1), s_max+1+s), 1), flip(prob(1:(s_max+1), s_max+1-s), 1)];
		
		enter_Es = [probE(:, s_max+1+s), probE(:, s_max+1-s)];

		enter_exogs = [prob_exog(:, s_max+1+s), prob_exog(:, s_max+1-s)];

		enter_Ls(s+1, :) = 0;
		enter_Fs(s+1, :) = 0;
		enter_Es(s+1, :) = 0;
		enter_exogs(s, :) = 0;

		H_Upsilon(rownum, 3:(s_max+2)) = flip(firm_distributions(2:end)).*flip(enter_Fs(2:end, 1) + enter_Fs(2:end, 2), 1)';
		H_Upsilon(rownum, s_max+3) = firm_distributions(1)*(enter_Fs(1, 1) + enter_Fs(1,2) + enter_Ls(1,1) + enter_Ls(1,2));
		H_Upsilon(rownum, (s_max+4):(2*s_max+3)) = firm_distributions(2:end).*(enter_Ls(2:end, 1) + enter_Ls(2:end, 2))';

		H_Upsilon(rownum, 2*s_max+4) = -(1-probE(s+1, s_max+1+s))*firm_distributions(s+1) + ...
										sum(firm_distributions.*(enter_Es(:,1) + enter_Es(:,2))'); 


		H_Upsilon(rownum, s_max+3+s) = H_Upsilon(rownum, s_max+3+s) - firm_distributions(s+1);
		H_Upsilon(rownum, s_max+3-s) = H_Upsilon(rownum, s_max+3-s) - firm_distributions(s+1)*(1-prob(s_max+1-s, s_max+1+s));

		H_Upsilon(rownum, 2*s_max+5) = investment(s_max+1)*(enter_Ls(1,1) + enter_Ls(1,2) + enter_Fs(1,1) + enter_Fs(1,2)) + ...
									   investment_entrant(1)*(enter_Es(1,1) + enter_Es(1,2));




		for (ii = 1:s_max)
			H_Upsilon(rownum, 2*s_max+5+ii) = investment(s_max+1-ii)*(enter_Fs(ii+1, 1) + enter_Fs(ii+1, 2)) + ...
											  investment(s_max+1+ii)*(enter_Ls(ii+1	,1) + enter_Ls(ii+1,2)) + ...
											  investment_entrant(1+ii)*(enter_Es(ii+1,1) + enter_Es(ii+1,2)) + ...
											  eta(ii+1)*(enter_exogs(ii, 1) + enter_exogs(ii, 2));
		end

		H_Upsilon(rownum, 2*s_max+5+s) = -((investment(s_max+1+s) + investment(s_max+1-s)*(1-prob(s_max+1-s, s_max+1+s)) + ...
			investment_entrant(s+1)*(1-probE(s+1, s_max+1+s)) + ...
			eta(s+1))); 
	end
	H_Upsilon(s_max+3, (2*s_max+5):end) = 1;

	column_index = -s_max:s_max;

	constant_func = @(B, gamma) B^(gamma/(1-gamma)) * gamma^(gamma/(1-gamma)) * B;
	for (ii = 1:length(-s_max:s_max))

		rownum = s_max+3+ii;
		s = column_index(ii);

		H_Upsilon(rownum, 2) = ...
			...% cs(investment(ii), B, subsidy(ii))*
			(sum(prob(ii, (ii+1):end).*(parv_paromega((ii+1):end) - parv_paromega(ii)) ) - ...
													Dv_omega(ii));
	end

	for (rownum = (s_max+4):(3*s_max+4))
		s = rownum - (2*s_max+4);
		idx = find(s == column_index);
		for (sigma = 1:length(-s_max:s_max))
			parv_parxsig_iter = parv_parxsig(sigma, :);
			H_Upsilon(rownum, 2+sigma) = ...
				 cs(investment(idx), B, subsidy(idx))*sum(prob(idx, (idx+1):end).*(parv_parxsig_iter((idx+1):end) - ...
					 																   parv_parxsig_iter(idx)));

		end


		H_Upsilon(rownum, 2*s_max+4) = cs(investment(idx), B, subsidy(idx))*sum(prob(idx, (idx+1):end).*...
											(parv_parxe(1, (idx+1):end) - parv_parxe(1, idx)));



	end


	for (ii = 1:length(-s_max:s_max))
		rownum = s_max+3+ii;
		s = column_index(ii);

		H_Upsilon(rownum, ii+2) = H_Upsilon(rownum, ii+2) - 1;
	end

	
	if (strcmp(entry_type, "Undirected"))

		other_term = zeros(1, s_max+1);
		other_termv = zeros(1, s_max+1);
		for (sigma = 1:length(other_term))
			other_term(sigma) = sum(probE(sigma, (s_max+3-sigma):end).*parv_paromega((s_max+3-sigma):end));
			other_termv(sigma) = sum(probE(sigma, (s_max+3-sigma):end).*value_functions((s_max+3-sigma):end));
		end

		H_Upsilon(3*s_max+5, 2) = cs(investment_entrant(1), B_entrant, subsidyE(1))*(sum(firm_distributions.*other_term) - ...
			sum(firm_distributions.*other_termv)/wage_share); 

		for (sigma = 1:length(-s_max:s_max))
			parv_parxsig_iter = parv_parxsig(sigma, :);

			other_term = zeros(1, s_max+1);
			for (jj = 1:length(0:s_max))
				other_term(jj) = sum(probE(jj, (s_max+3-jj):end).*parv_parxsig_iter((s_max+3-jj):end));
			end

			H_Upsilon(3*s_max+5, 2+sigma) = ...
				cs(investment_entrant(1), B_entrant, subsidyE(1))*sum(firm_distributions.*other_term);

		end

		other_term = zeros(1, s_max+1);
		for (jj = 1:length(other_term))
			other_term(jj) = sum(probE(jj, (s_max+3-jj):end).*parv_parxe(1, (s_max+3-jj):end));
		end

		H_Upsilon(3*s_max+5, 2*s_max+4) = ...
			cs(investment_entrant(1), B_entrant, subsidyE(1))*sum(firm_distributions.*other_term) - 1;


		for (sigma = 1:length(0:s_max))
  			H_Upsilon(3*s_max+5, 2*s_max+4+sigma) = ...
				cs(investment_entrant(1), B_entrant, subsidyE(1))*sum(probE(sigma, (s_max+3-sigma):end).*value_functions((s_max+3-sigma):end));
		end
			
	end

	row_remove = [];
	column_remove = [];
	if (elastic_labor)
		row_remove = [row_remove, 1];
		column_remove = [column_remove, 2];
	end

	if (~entrance) 
		row_remove = [row_remove, 3*s_max+5];
		column_remove = [column_remove, 2*s_max+4];
	end

	H_Upsilon(row_remove, :) = [];
	H_Upsilon(:, column_remove) = [];


	divisor_I = eye(2*s_max+1 + 1*(entrance));
	divisor_zero = zeros(size(H_Upsilon,1) - size(divisor_I,2),size(divisor_I,2));
	divisor = [divisor_zero; divisor_I];

	if (abs(m.gamma - 0.5) > 1e-6)
		[r,c] = find(isnan(H_Upsilon)); 
		H_Upsilon(unique(r), :) = [];

		divisor(unique(r), :) = [];
	end

	M = -H_Upsilon \ divisor;
	if (~entrance) && (elastic_labor)
		if (strcmp(block, "all"))
			M = -H_Upsilon \ divisor;
		elseif (strcmp(block, "all_but_g"))
			H_iter = H_Upsilon;
			H_iter(1,:) = [];
			H_iter(:,1) = [];

			divisor_iter = [zeros(size(H_iter,1) - size(divisor_I,2),size(divisor_I,2)); divisor_I];
			M = -H_iter \ divisor_iter;
		elseif (strcmp(block, "innovation"))
			H_iter = H_Upsilon((s_max+3):(3*s_max+3), 2:(2*s_max+2));
			divisor_iter = [zeros(size(H_iter,1) - size(divisor_I,2),size(divisor_I,2)); divisor_I];
			M = -H_iter \ divisor_iter;
		end
	end
	m.M = M; 

end