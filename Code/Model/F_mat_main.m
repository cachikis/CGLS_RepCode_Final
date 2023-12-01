%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: F_mat_main.m
% Author: Craig A. Chikis
% Date: 11/30/2020
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prob, probE, prob_exog] = F_mat_main(phi_wt, phi_wt_tilde, phi_wt_exog, l, l_tilde, s_max)





	prob = zeros(2*s_max+1);
	probE = zeros(s_max+1, 2*s_max+1);
	prob_exog = zeros(s_max, 2*s_max+1);

	for (ii = 1:(s_max-1))
		prob(s_max+1+ii, max(s_max+1+ii+1, s_max+1+l)) = prob(s_max+1+ii, max(s_max+1+ii+1, s_max+1+l)) + phi_wt;
		prob(s_max+1+ii, s_max+1+ii+1) = prob(s_max+1+ii, s_max+1+ii+1) + 1-phi_wt;

		
		probE(ii+1, s_max+1+l_tilde) = probE(ii+1, s_max+1+l_tilde) + phi_wt_tilde;
		probE(ii+1, s_max+1-ii+1) = probE(ii+1, s_max+1-ii+1) + 1-phi_wt_tilde;

		prob_exog(ii, s_max+1) = prob_exog(ii, s_max+1) + phi_wt_exog;
		prob_exog(ii, s_max+1-ii+1) = prob_exog(ii, s_max+1-ii+1) + 1-phi_wt_exog;

		prob(s_max+1-ii, s_max+1+l) = prob(s_max+1-ii, s_max+1+l) + phi_wt;
		prob(s_max+1-ii, s_max+1-ii+1) = prob(s_max+1-ii, s_max+1-ii+1) + 1-phi_wt;
	end
	prob(s_max+1, max(s_max+1+1, s_max+1+l)) = prob(s_max+1, max(s_max+1+1, s_max+1+l)) + phi_wt;
	prob(s_max+1, s_max+1+1) = prob(s_max+1, s_max+1+1) + 1-phi_wt;

	prob(2*s_max+1,2*s_max+1) = 1;

	prob(1, 2) = prob(1,2) + 1-phi_wt;
	prob(1, s_max+1+l) = prob(1, s_max+1) + phi_wt;

	probE(1, max(s_max+1+l_tilde, s_max+1+1)) = probE(1, max(s_max+1+l_tilde, s_max+1+1)) + phi_wt_tilde;
	probE(1, s_max+1+1) = probE(1, s_max+1+1) + 1-phi_wt_tilde;

	probE(s_max+1, s_max+1+l_tilde) = probE(s_max+1, s_max+1+l_tilde) + phi_wt_tilde;
	probE(s_max+1, 2) = probE(s_max+1, 2) + 1-phi_wt_tilde;

	prob_exog(s_max, s_max+1) = prob_exog(s_max, s_max+1) + phi_wt_exog;
	prob_exog(s_max, 2) = prob_exog(s_max, 2) + 1-phi_wt_exog;

	prob_exog = flip(prob_exog,1);
	probE = flip(probE,1);



end
