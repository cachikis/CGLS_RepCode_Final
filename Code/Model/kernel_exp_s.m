%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: kernel_exp_s.m
% Author: Craig A. Chikis
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prob, probE, prob_exog] = kernel_exp_s(phi_bar, zeta_bar, phi_E_bar,...
                                                 nu_phi, nu_zeta, nu_phi_E, ...
                                                 s_max, jump_max, jump_max_p, jump_max_E)

    prob = zeros(2*s_max+1);
    probE = zeros(s_max+1, 2*s_max+1);
    prob_exog = zeros(s_max, 2*s_max+1);

    phi_func = @(phi_bar, s, parameter) phi_bar*exp(parameter*(1-s));

    for (ii = 1:(s_max-1))
       
        prob(ii, ii+1) = 1 - phi_func(phi_bar, s_max+1-ii, nu_phi);
        prob(ii, ii+jump_max) = prob(ii, ii+jump_max) + phi_func(phi_bar, s_max+1-ii, nu_phi);

        prob(s_max+ii, s_max+ii+1) = 1;
    end
    prob(1:s_max, s_max+1) = sum(prob(1:s_max, (s_max+2):end), 2) + prob(1:s_max, s_max+1);
    prob(1:s_max, (s_max+2):end) = 0;

    prob(s_max, s_max+1) = 1;
    prob(2*s_max+1, 2*s_max+1) = 1;
    prob(2*s_max, 2*s_max+1) = 1;

    prob = prob./sum(prob, 2);

    for (ii = 1:size(prob_exog, 1))
        prob_exog(ii, ii+1) = 1 - phi_func(zeta_bar, s_max+1-ii, nu_zeta);
        prob_exog(ii, ii+jump_max_p) = prob_exog(ii, ii+jump_max_p) + phi_func(zeta_bar, s_max+1-ii, nu_zeta);
    end
    prob_exog(1:s_max, s_max+1) = sum(prob_exog(1:s_max, (s_max+2):end), 2) + prob_exog(1:s_max, s_max+1);
    prob_exog(1:s_max, (s_max+2):end) = 0;

    for (ii = 1:(size(probE, 1)-1))
        probE(ii, ii+1) = 1 - phi_func(phi_E_bar, s_max+2-ii, nu_phi_E);
        probE(ii, ii+jump_max_E) = probE(ii, ii+jump_max_E) + phi_func(phi_E_bar, s_max+2-ii, nu_phi_E);
    end
    probE(1:s_max, s_max+1) = sum(probE(1:s_max, (s_max+2):end), 2) + probE(1:s_max, s_max+1);
    probE(1:s_max, (s_max+2):end) = 0;
    probE(s_max+1, (s_max+2):end) = 0;
    probE(s_max+1, s_max+2) = 1;

    prob_exog = prob_exog./sum(prob_exog, 2);
    probE=  probE./sum(probE, 2);

end
