%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: firm_distributions_v3.m
% Author: Craig A. Chikis
% Date: 06/08/2020
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, mu, b] = firm_distributions_v3(investment, investmentE, eta, probL, probF, prob0, ...
 										 probE, prob_exog, s_max)
	

	


	xL = investment((s_max+2):end);
	xF = flip(investment(1:s_max));
	x0 = investment(s_max+1);

	xE = investmentE;

	A = zeros(s_max+1, s_max+1);
	b = zeros(s_max+1, 1);

	b(end) = 1;
	A(s_max+1, :) = 1;

	prob = [flip(probF, 1);prob0;probL];
	
	prob0 = prob(:, s_max+1);
	prob0_F = flip(prob0(1:s_max))';

	prob0_exog = prob_exog(:, s_max+1);
	prob0_exog_F = prob0_exog(1:s_max)';

	prob0_E = probE(2:end, s_max+1)';
	prob0_E_F = prob0_E(1:s_max);

	dot_F = (prob0_F.*xF + prob0_exog_F.*eta(2:end) + prob0_E_F.*xE(2:end));

	A(1, 2:end) = A(1, 2:end) + dot_F;

	A(1, 1) = -(2*x0 +xE(1)*(1-probE(1, s_max+1))); 

	for (i = 2:(size(A,1)-1))
		A(i, i) = -(xL(i-1) + xF(i-1)*(1 - prob(s_max+1-(i-1), s_max+1+(i-1))) + eta(i) + xE(i)*(1 - probE(i, s_max+i))); 
	end


	for (i = 2:(size(A,1)-1))
		prob_enter = [prob(:, (s_max+1)+i-1), prob(:, (s_max+1)-i+1)];

		prob_enter_F = flip(prob_enter(1:s_max, :), 1);
		prob_enter_F(i-1, :) = 0;

		prob_enter_L = prob_enter((s_max+2):end, :);
		prob_enter_L(i-1, :) = 0;

		prob_enter_0 = prob_enter(s_max+1, :);


		probE_enter = [probE(:, (s_max+1)+i-1), probE(:, (s_max+1)-i+1)];
		probE_enter_0 = probE_enter(1, :);

		probE_enter_F = probE_enter(2:end, :);
		probE_enter_F(i-1, :) = 0;

		prob_exog_enter = [prob_exog(:, (s_max+1)+i-1), prob_exog(:, (s_max+1)-i+1)];

		dot_F = (prob_enter_L(:, 1) + prob_enter_L(:, 2))'.*xL + ...
				(prob_enter_F(:, 1) + prob_enter_F(:, 2))'.*xF + ...
				(probE_enter_F(:, 1) + probE_enter_F(:, 2))'.*xE(2:end) + ...
				(prob_exog_enter(:, 1) + prob_exog_enter(:,2))'.*eta(2:end);

		A(i, 2:end) = A(i, 2:end) + dot_F;
		A(i, 1) = 2*x0*(prob_enter_0(:, 1) + prob_enter_0(:, 2)) + ...
					xE(1)*(probE_enter_0(:,1) + probE_enter_0(:,2));

	end




	mu = A\b;


end