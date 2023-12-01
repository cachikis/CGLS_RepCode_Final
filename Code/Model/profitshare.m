%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: profitshare.m
% Author: Craig A. Chikis
% Date: 09/03/2023
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = profitshare(m)

    prof = m.profit ./ (1 - m.tax_rate); 

    m.profitshare = sum(m.firm_distributions .* (m.profit((m.s_max+1):end) + flip(m.profit(1:(m.s_max+1))))); 


    
end