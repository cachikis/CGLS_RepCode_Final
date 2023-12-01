%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: F_mat_AAI.m
% Author: Craig A. Chikis
% Date: 07/29/2022
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prob, probE, prob_exog] = F_mat_AAI(phi_wt, phi_wt_tilde, phi_wt_exog, l, l_tilde, s_max)


    F = @(phi_wt, s, s_max) (s + s_max).^(-(-phi_wt));

	prob = zeros(2*s_max+1);
	probE = zeros(s_max+1, 2*s_max+1);
	prob_exog = zeros(s_max, 2*s_max+1);


    prob(2*s_max+1, 2*s_max+1) = 1;
    srange = -s_max:s_max;
    prob(1, 2:(s_max+1+l)) = prob(1, 2:(s_max+1+l)) + F(phi_wt, (srange(1)+1):l, s_max);
    prob(1, :) = prob(1, :)/sum(prob(1, :));
    for (ii = 2:(2*s_max))
        if (srange(ii) + 1 < l)
            prob(ii, ii+1) = F(phi_wt, srange(ii) + 1, s_max) + sum(F(phi_wt , (-s_max + 1):(srange(ii)), s_max));
            prob(ii, (ii+2):(s_max+1+l)) = F(phi_wt, (srange(ii)+2):l, s_max);
        else
            prob(ii, ii+1) = 1;
        end
    end
    prob(2:end, :) = prob(2:end, :)./sum(prob(2:end, :), 2);

	for (ii = 1:(s_max-1))
	

		
		probE(ii+1, s_max+1+l_tilde) = probE(ii+1, s_max+1+l_tilde) + phi_wt_tilde;
		probE(ii+1, s_max+1-ii+1) = probE(ii+1, s_max+1-ii+1) + 1-phi_wt_tilde;

		prob_exog(ii, s_max+1) = prob_exog(ii, s_max+1) + phi_wt_exog;
		prob_exog(ii, s_max+1-ii+1) = prob_exog(ii, s_max+1-ii+1) + 1-phi_wt_exog;

	
	end


	probE(1, max(s_max+1+l_tilde, s_max+1+1)) = probE(1, max(s_max+1+l_tilde, s_max+1+1)) + phi_wt_tilde;
	probE(1, s_max+1+1) = probE(1, s_max+1+1) + 1-phi_wt_tilde;

	probE(s_max+1, s_max+1+l_tilde) = probE(s_max+1, s_max+1+l_tilde) + phi_wt_tilde;
	probE(s_max+1, 2) = probE(s_max+1, 2) + 1-phi_wt_tilde;

	prob_exog(s_max, s_max+1) = prob_exog(s_max, s_max+1) + phi_wt_exog;
	prob_exog(s_max, 2) = prob_exog(s_max, 2) + 1-phi_wt_exog;

	prob_exog = flip(prob_exog,1);
	probE = flip(probE,1);



end
