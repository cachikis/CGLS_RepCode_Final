%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: hhi.m
% Author: Craig A. Chikis
% Date: 09/03/2023
% Note(s): 
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = hhi(m)

    % sales_divisor = sales_func(m.kappa, m.s_max, m.nu_s); 
    % m.panel_save(m.panel_save.Industry == 1, :)
    tmp1 = groupsummary(m.panel_save(:, ["Firm", "Industry", "Time", "Sale_t"]), ["Firm", "Industry", "Time"], 'sum');
    tmp2 = groupsummary(m.panel_save(:, ["Industry", "Time", "Sale_t"]), ["Industry", "Time"], 'sum');
    tmp1 = innerjoin(tmp1, tmp2, 'keys', ["Industry", "Time"]);
    tmp1.share = tmp1.sum_Sale_t_tmp1 ./ tmp1.sum_Sale_t_tmp2; 
    tmp3 = groupsummary(tmp1(:, ["Industry", "Time", "share"]), ["Industry", "Time"], @(x) sum(x .^ 2));

    m.hhi = mean(tmp3.fun1_share);
    % m.hhi = sum(m.firm_distributions .* (flip(sales_divisor(1:(m.s_max+1))).^2 + sales_divisor((m.s_max+1):end).^2));


end