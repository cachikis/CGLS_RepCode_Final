%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: transition_constant_surprise2.m
% Author: Craig A. Chikis
% Date: 11/19/2023
% Note(s)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x_forward, v_forward] = transition_constant_surprise(m, rhostart, rhoend, ...
                                                               years, yearsrhochange, ...
                                                               justprodlabor, feedmu, ...
                                                               rhotilde, rho_vec, ...
                                                               rho_t)

    s_max = m.s_max; 

    if (~m.elastic_labor) 
        error('Prod. labor not built for inelastic labor.')
    end

    G = @(x,B,gamma) (x./B).^(1./gamma);
    GP = @(x,B,gamma) (1./(B.*gamma)).*(x./B).^(1./gamma - 1);
    GPInv = @(z,B,gamma) B.*(z.*B.*gamma).^(gamma./(1-gamma));
                                

                                
    % rho_vec = power(1 + [rhostart, rhoend], 1/12) - 1;
    % rho_t = [repelem(rhotilde+power(1+rhostart, 1/12)-1, 12*years), ...
    %         power(1 + linspace(rhostart+(1+rhotilde)^12-1+1e-6, rhoend+(1+rhotilde)^12-1, 12*yearsrhochange), 1/12)-1, ...
    %         repelem(power(1+rhoend+(1+rhotilde)^12-1,1/12)-1, 12*years)];

    v_collect = zeros(length(rho_vec), 2*s_max+1);
    mu_collect = zeros(length(rho_vec), s_max+1);
    x_collect = zeros(size(v_collect));
    profit_collect = zeros(size(v_collect));

    for (ii = 1:length(rho_vec))

        m_tmp = m; 
        m_tmp.rho = rho_vec(ii); %rhotilde + rho_vec(ii); 
        m_tmp = workhorse_nested_robust(m_tmp); 

        v_collect(ii, :) = m_tmp.value_functions;
        x_collect(ii, :) = m_tmp.investment;
        mu_collect(ii, :) = m_tmp.firm_distributions;
        profit_collect(ii, :) = m_tmp.profit; 

    end

    if (~all(profit_collect(1, :) == profit_collect(2, :))) 
        error('Need to code this for time variation in profit.')
    end
    if (m.entrance) 
        error('See below for shoot forward. Code not set up for shooting forward.')
    end

    rho_perm = rho_t((12*years+1):(12*years+1+yearsrhochange-1));
    x_trans = zeros(length(rho_perm), 2*s_max+1);
    v_trans = zeros(size(x_trans));
    m_cell = cell(size(rho_perm));  
    for (ii = 1:length(rho_perm))
        m_tmp = m; 
        m_tmp.rho = rho_perm(ii); 
        m_cell{ii} = m_tmp; 
    end
    parfor (ii = 1:length(rho_perm))

            m_cell{ii} = workhorse_nested_robust(m_cell{ii});

            v_trans(ii, :) = m_cell{ii}.value_functions;
            x_trans(ii, :) = m_cell{ii}.investment;

    end
    x_forward = [repmat(x_collect(1, :), [12*years, 1]); x_trans; repmat(x_collect(2, :), [12*years, 1]); nan(1,2*s_max+1)]; 
    v_forward = [repmat(v_collect(1, :), [12*years, 1]); v_trans; repmat(v_collect(2, :), [12*years, 1]); v_collect(2,:)]; 





end