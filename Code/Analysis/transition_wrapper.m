%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: transition_wrapper.m
% Author: Craig A. Chikis
% Date: 09/12/2022
% Note(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [prod, quality, dL_t_dt, ...
          muder, xder, rdlabor, prodlabor, ...
          prodlaborder, markup_vec, ...
          rho_path, tplot, ximpact, ...
          muimpact, ximpactback, muimpactback, tplotback, ...
          quality2, ...
          prod2, ...
          prod3] = transition_wrapper(m, rhostart, rhoend, years, yearsrhochange, justprodlabor, feedmu, no_surprise, rhotilde) 

    s_max = m.s_max; 
    prob = m.prob; 
    eta = m.eta; 
    prob_exog = m.prob_exog;
    probE = m.probE; 
    kappa = m.kappa; 
    subsidy = m.subsidy; 

    
    if (~m.elastic_labor)
        error('Prod. labor not built for this. Also \omega \ne 1 so FOCs for GPInv wrong, among other things.')
    end

    if (any(sum(prob(1:s_max, (s_max+2):end), 2) > 0)) || (abs(prob(s_max+1, s_max+2) - 1) >= 1e-8) 
        error("quality index growth can't handle leapfrogging")
    end

    G = @(x,B,gamma,subsidy) (1-subsidy) .* (x./B).^(1./gamma);
    GP = @(x,B,gamma,subsidy) (1-subsidy) .* (1./(B.*gamma)).*(x./B).^(1./gamma - 1);
    GPInv = @(z,B,gamma,subsidy) B.*(z.*B.*gamma ./ (1-subsidy)).^(gamma./(1-gamma));
                                

    rho_vec = power(1 + [rhostart, rhoend], 1/12) - 1;
    rho_t = [repelem(rhotilde+power(1+rhostart, 1/12)-1, 12*years), ...
            linspace(rhotilde+(power(1+rhostart, 1/12)-1)+1e-10, rhotilde+(power(1+rhoend, 1/12)-1), 12*yearsrhochange), ...
            ...% power(1 + linspace(rhostart+(1+rhotilde)^12-1+1e-6, rhoend+(1+rhotilde)^12-1, 12*yearsrhochange), 1/12)-1, ...
            repelem(rhotilde+power(1+rhoend, 1/12)-1, 12*years)];

    
    v_collect = zeros(length(rho_vec), 2*s_max+1);
    mu_collect = zeros(length(rho_vec), s_max+1);
    x_collect = zeros(size(v_collect));
    profit_collect = zeros(size(v_collect));
    lab_collect = cell(1, length(rho_vec));

    for (ii = 1:length(rho_vec))
        m_tmp = m; 
        m_tmp.rho = rhotilde + rho_vec(ii); 
        m_tmp = workhorse_nested_robust(m_tmp); 


        v_collect(ii, :) = m_tmp.value_functions;
        x_collect(ii, :) = m_tmp.investment;
        mu_collect(ii, :) = m_tmp.firm_distributions;
        profit_collect(ii, :) = m_tmp.profit; 
        lab_collect{ii} = m_tmp.lab_demand; 

    end
    B = m_tmp.B;
    gamma = m_tmp.gamma; 
    lambda = m_tmp.lambda;
    nu_s = m_tmp.nu_s; 

    if ( (~all(profit_collect(1, :) == profit_collect(2, :))) || (~all(lab_collect{1}{2} == lab_collect{2}{2})) )
        error('Need to code this for time variation in profit.')
    end
    if (m.entrance)
        error('See below for shoot forward. Code not set up for shooting forward with entry or inelastic labor.')
    end

    x_forward = zeros(length(rho_t)+1, 2*s_max+1);
    v_forward = zeros(length(rho_t), 2*s_max+1); 
    del_t = 1;
    rhoflip = flip(rho_t);

    v_forward(1, :) = v_collect(2, :); 
    x_forward(1, :) = nan;
    x_forward(2, :) = GPInv(sum(prob.*v_forward(1, :), 2)' - v_forward(1, :), B, gamma, subsidy);

    prob_rearrange = [flip(prob((s_max+2):end, :), 1); prob(s_max+1, :); flip(prob(1:s_max, :), 1)];
    prob_rearrange = flip(prob_rearrange, 2);

    prob_exog_rearrange = [flip(prob_exog(1:s_max, :), 1)];
    prob_exog_rearrange = flip(prob_exog_rearrange, 2);
    if (no_surprise)
        for (ii = 2:length(rho_t))

            x_c = [flip(x_forward(ii, (s_max+2):end)), x_forward(ii, s_max+1), flip(x_forward(ii, 1:s_max))];

            RHS = profit_collect(1, :) - ...
                        G(x_forward(ii, :), B, gamma,subsidy) + ...
                        x_forward(ii, :).*(sum(prob.*v_forward(ii-1, :), 2)' - v_forward(ii-1, :)) + ...
                        x_c.*(sum(prob_rearrange.*v_forward(ii-1, :), 2)' - v_forward(ii-1, :));

            RHS(1:s_max) = RHS(1:s_max) + ...
                flip(eta(2:end)).*(sum(prob_exog.*v_forward(ii-1, :), 2)' - v_forward(ii-1, 1:s_max));

            RHS((s_max+2):end) = RHS((s_max+2):end) + ...
                eta(2:end).*(sum(prob_exog_rearrange.*v_forward(ii-1, :), 2)' - v_forward(ii-1, (s_max+2):end)); 
            


            v_forward(ii, :) = (RHS + v_forward(ii-1, :)/del_t)./(rhoflip(ii) + 1/del_t); 


            x_tmp = max(GPInv(sum(prob.*v_forward(ii, :), 2)' - v_forward(ii, :), B, gamma, subsidy), 0);
            x_tmp(imag(x_tmp) ~= 0) = 0;
            x_forward(ii+1, :) = x_tmp; 

        end
        x_forward = flip(x_forward, 1);
        v_forward = flip(v_forward, 1);
    else
         [x_forward, v_forward] = transition_constant_surprise(m, rhostart, rhoend, years, yearsrhochange, justprodlabor, feedmu, rhotilde);
    end
  
    mu_forward = [nan(12*years, s_max+1); zeros(length(rho_t)-12*years, s_max+1)];
    if (length(feedmu) ~= s_max+1)
        mu_forward(12*years+1, :) = mu_collect(1, :); 
    else
        mu_forward(12*years+1, :) = feedmu;
    end

    for (ii = 1:(size(mu_forward((12*years+1):end, :), 1)-1))
        M = transmat_matrix_instantiation(prob, probE, prob_exog, x_forward(12*years+ii, :), zeros(1, s_max+1), eta, s_max);
        mu_forward(12*years+ii+1, :) = mu_forward(12*years+ii, :)*M; 
    end


    rho_path = 100*(power(1 + [rhotilde + power(1+rhostart,1/12)-1, rho_t((12*years+1):end)], 12) - 1);
    
    tplot = (0:(12*yearsrhochange + 12*years))/12; 
    
    ximpact = [x_collect(1, :);
               x_forward((12*years+1):(end-1), :)]; 
    
    muimpact = [mu_forward(12*years+1, :);
                mu_forward((12*years+2):end, :);
                nan(1, size(mu_forward, 2))];
    
    ximpactback = [x_forward(1:(12*years + 12*yearsrhochange), :)];
    
    muimpactback = [mu_forward(1:(12*years + 12*yearsrhochange), :)];
    
    tplotback = [(-12*years:-1)/12, (0:(12*yearsrhochange-1))/12];


    quality = zeros(1, size(ximpact, 1));
    quality2 = zeros(size(quality)); 
    prod = zeros(size(quality));
    prod2 = nan(size(quality)); 
    prod3 = nan(size(quality));
    term1 = nan(size(quality));
    term2 = nan(size(quality)); 
    muder = zeros(size(ximpact, 1), s_max+1);
    xder = zeros(size(ximpact, 1), 2*s_max+1);
    rdlabor = zeros(size(xder));
    prodlabor = zeros(size(xder));
    prodlaborder = zeros(size(xder));
    markup_vec = zeros(3, size(quality, 2));

    for (ii = 1:size(ximpact, 1))

        try 
            muder(ii, :) = (muimpact(ii, :) - muimpact(ii-1, :))/(del_t);
            xder(ii, :) = (ximpact(ii, :) - ximpact(ii-1, :))/(del_t);
        catch
            muder(ii, :) = nan;
            xder(ii, :) = nan;
        end

        rdlabor(ii, :) = G(ximpact(ii, :), B, gamma, 0); 
        prodlabor(ii, :) = lab_collect{1}{2}; 
        try
            [p90_markup, markup_array, ~] =  markup_quant(kappa, lambda, s_max, muimpact(ii, :), nu_s, 0.9);
            [p50_markup, markup_array, ~] =  markup_quant(kappa, lambda, s_max, muimpact(ii, :), nu_s, 0.5);
            markup_vec(1, ii) = 100*(sum(markup_array(1, :).*muimpact(ii, :)) - 1);
            markup_vec(2, ii) = 100*(p90_markup-1);
            markup_vec(3, ii) = 100*(p50_markup-1);
        catch 
            markup_vec(1, ii) = nan;
            markup_vec(2, ii) = nan;
            markup_vec(3, ii) = nan;
        end

    end
    mu0der = muder(:, 1); 
    
    dL_t_dt = zeros(1, size(xder, 1));
    for (ii = 1:size(xder, 1))

        try 
            prodlaborder(ii, :) = (prodlabor(ii, :) - prodlabor(ii-1, :))/(del_t);
        catch
            prodlaborder(ii, :) = nan;
        end

        L_t_nosum = (~justprodlabor)*(rdlabor(ii, (s_max+1):end) + flip(rdlabor(ii, 1:(s_max+1)))) + ...
                                        prodlabor(ii, (s_max+1):end) + flip(prodlabor(ii, 1:(s_max+1)));

        L_t = sum(muimpact(ii, :).*L_t_nosum);  
        dL_t_dt(ii) = (1./L_t).*( sum(muder(ii, :).*L_t_nosum) + ...
                                (~justprodlabor)*sum(muimpact(ii, :).*(GP(ximpact(ii, (s_max+1):end), B, gamma, subsidy((s_max+1):end)).*xder(ii, (s_max+1):end) + ...
                                                    GP(flip(ximpact(ii, 1:(s_max+1))), B, gamma, flip(subsidy(1:(s_max+1)))).*flip(xder(ii, 1:(s_max+1))) + ...
                                                    prodlaborder(ii, (s_max+1):end) + flip(prodlaborder(ii, 1:(s_max+1))) )) ); 


        quality(ii) = log(lambda)*sum(muimpact(ii, :).*ximpact(ii, (s_max+1):end).*[2, ones(1, s_max)]); 
        quality2(ii) = log(lambda)*sum(muimpact(ii, 2:end).*(flip(ximpact(ii, 1:s_max)) + eta(2:end))); 
        prod(ii) = quality(ii) - log(lambda)*sum((0:s_max).*muder(ii, :)) - dL_t_dt(ii);

        term1(ii) =  mu0der(ii) * ( (kappa/(kappa-1))*log((lambda * prodlabor(ii, s_max+2))^((kappa-1)/kappa) + ...
        (prodlabor(ii, s_max))^((kappa-1)/kappa) ) - ...
(kappa/(kappa-1))*log(2 * prodlabor(ii, s_max+1)^((kappa-1)/kappa) ) ); 

        term2(ii) =   mu0der(ii) * ((prodlabor(ii, s_max+2) + prodlabor(ii, s_max) - 2*prodlabor(ii, s_max+1))/( ...
        (1-muimpact(ii,1))*(prodlabor(ii, s_max+2)+prodlabor(ii,s_max)) + ...
            2*muimpact(ii,1)*prodlabor(ii,s_max+1) ) );  

        prod2(ii) = quality2(ii) - term1(ii) + term2(ii); 
        prod3(ii) = quality(ii) - term1(ii) + term2(ii); 
                  
    end 
    prod(1) = log(lambda)*sum(ximpact(1, (s_max+1):end).*[2, ones(1,s_max)].*muimpact(1, :));
    prod2(1) = prod(1);
    prod3(1) = prod(1); 
    

end