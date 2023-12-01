%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: transmat_matrix_instantiation.m
% Author: Craig A. Chikis
% Date: 09/10/2020
% MA-MFA
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function M = transmat_matrix_instantiation(prob, probE, prob_exog, investment, investment_entrant, eta, s_max)


	M = zeros(s_max+1);
	prob_exog = [prob_exog; zeros(1, size(prob_exog,2))];

	for (ii = 0:s_max)

			probL_iter = prob(s_max+1+ii, :);

			probF_iter = prob(s_max+1-ii, :);
			probF_iter(s_max+1+ii) = 0;

			prob_exog_iter = prob_exog(s_max+1-ii, :);

			probE_iter = probE(s_max+1-ii, :);
			probE_iter(s_max+1+ii) = 0;


		
			leaders_ahead = probL_iter((s_max+2):end)*investment(s_max+1+ii);

			followers_no_leap = flip(probF_iter(1:s_max))*investment(s_max+1-ii);
			followers_leap = probF_iter((s_max+2):end)*investment(s_max+1-ii);
			followers_tied = probF_iter(s_max+1)*investment(s_max+1-ii);

			exog_no_leap = flip(prob_exog_iter(1:s_max))*eta(ii+1);
			exog_tied = prob_exog_iter(s_max+1)*eta(ii+1);

			entrant_no_leap = flip(probE_iter(1:s_max))*investment_entrant(ii+1);
			entrant_leap = probE_iter((s_max+2):end)*investment_entrant(ii+1);
			entrant_tied = probE_iter(s_max+1)*investment_entrant(ii+1);


			if (ii > 0)
				followers_leap(ii) = 0;
				entrant_leap(ii) = 0;
			end


			M(ii+1, 2:end) = leaders_ahead + followers_no_leap + followers_leap + exog_no_leap + ...
						     entrant_no_leap + entrant_leap;
			M(ii+1, 1) = followers_tied + exog_tied + entrant_tied;

			M(ii+1, ii+1) = M(ii+1, ii+1) + (1 - sum(M(ii+1, :)));

	end


end