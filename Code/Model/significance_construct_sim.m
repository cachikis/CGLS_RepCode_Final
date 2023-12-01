%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: significance_construct_sim.m
% Author: Craig A. Chikis
% Date: 10/18/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = significance_construct_sim(m, method_prc)

    if nargin < 2
        method_prc = "exact";
    end

    innovation_mat = m.innovation_mat; 
    zeta = m.zeta;
    lambda = m.lambda; 
    fcitcount = m.fcitcount;
    month_increment = m.month_increment;
    numCPC = m.numCPC; 
    prctile_inp = m.prctile_inp; 
    l = m.l; 
    tech_cluster_prob = m.tech_cluster_prob;

    rng(20230903, 'Threefry');


    a = 1;

    stateL = innovation_mat(:, :, 1, 2);
    stateF = innovation_mat(:, :, 4, 2); 

    numfirms = 2*size(innovation_mat, 1) + sum(innovation_mat(:, :, [3,6,9,10], 1), 'all');
    total_patents = sum(sum(innovation_mat(:, :, :, 1), 3), 1);
    panel1 = zeros(numfirms/2, size(innovation_mat, 2)); 
    panel2 = zeros(size(panel1));
    panel3 = zeros(sum(innovation_mat(:, :, [3,9], 1), 'all'), size(innovation_mat, 2));
    panel4 = zeros(sum(innovation_mat(:, :, [6,10], 1), 'all'), size(innovation_mat, 2));


    [renterL, centerL] = find(sum(innovation_mat(:, :, [3,9], 1), 3) == 1);
    [renterF, centerF] = find(sum(innovation_mat(:, :, [6,10], 1), 3) == 1);
    [rleaderinc, cleaderinc] = find(sum(innovation_mat(:, :, [1,2,7], 1), 3) == 1);
    [rlaggardinc, claggardinc] = find(sum(innovation_mat(:, :, [4,5,8], 1), 3) == 1);

    panel1(sub2ind(size(panel1), rleaderinc, cleaderinc)) = 1;
    panel2(sub2ind(size(panel2), rlaggardinc, claggardinc)) = 1;
    panel3(sub2ind(size(panel3), renterL, centerL)) = 1;
    panel4(sub2ind(size(panel4), renterF, centerF)) = 1;


    panel = [panel1; panel2; panel3; panel4];
    statepanel = [stateL; stateF; stateL(1:size(panel3, 1), :); stateF(1:size(panel4, 1), :)];
    [rda, cda] = find(panel == 1);


    statepanel_reshape = reshape(statepanel', [month_increment, size(statepanel, 2)/month_increment, size(statepanel, 1)]);



    firmnum = [1:size(panel1, 1),  (size(panel1, 1)+1):(size(panel1, 1) + size(panel2, 1)), ...
                (size(panel1, 1) + size(panel2, 1)+1):(size(panel3, 1) + size(panel2, 1) + size(panel1, 1)), ...
                (size(panel1, 1) + size(panel2, 1) + size(panel3, 1)+1):(size(panel1, 1) + size(panel2, 1) + ...
                    size(panel3, 1) + size(panel4, 1))]'; 
    indpanel = [1:size(panel1, 1), 1:size(panel2, 1), renterL, renterF]';

    significance = zeros(length(rda), size(innovation_mat, 2));
    firmtrack = zeros(length(rda), 1);
    indtrack = zeros(length(rda), 1);
    changetrack = zeros(length(rda), 8);

    for (ii = 1:length(rda))

        firmtrack(ii) = firmnum(rda(ii));
        indtrack(ii) = indpanel(rda(ii));

        change = statepanel(rda(ii), cda(ii)) - ...
                    statepanel(rda(ii), cda(ii) - 1);
                

      

       
            if (~isnan(change))
                significance(ii, (cda(ii)):min(cda(ii)+fcitcount*month_increment-1, size(significance, 2))) = min(zeta*log(lambda^(change^a)), 1);
            else
                significance(ii, (cda(ii)):min(cda(ii)+fcitcount*month_increment-1, size(significance, 2)))  = nan;
            end

        changetrack(ii, :) = [firmtrack(ii), indtrack(ii), ...
                              statepanel(rda(ii), cda(ii)), statepanel(rda(ii), cda(ii)-1), ...
                              change, ii, cda(ii), max(statepanel_reshape(round(month_increment*(1 - (ceil(cda(ii)/month_increment) - cda(ii)/month_increment))):end, ...
                                                                           ceil(cda(ii)/month_increment), rda(ii)))];

    end
    

  
    uniquefirms = unique(firmtrack);
    CPCrand = randi(numCPC, length(unique(firmtrack)), 1);
    for (ii = 1:length(uniquefirms))
        indtrack(firmtrack == uniquefirms(ii)) = CPCrand(ii);
    end
    changetrack(:, 2) = indtrack;

    indtrack_unique = unique(indtrack); 
    for (ii = 1:length(indtrack_unique))
        changetrackiter = changetrack(changetrack(:, 2) == indtrack_unique(ii), :);

        [~, idx] = sort(changetrackiter(:, 7));
        changetrackiter = changetrackiter(idx, :); 

        newcluster = ...
                    ( ((changetrackiter(:, 4) < l) ) & ...
                        (changetrackiter(:, 3) - changetrackiter(:, 4) > 1) & ...
                            (changetrackiter(:, 3) >= l) & ...
                            (rand <= tech_cluster_prob)   ); 

          
        amend = changetrackiter(newcluster, :); 
   

        for (jj = 1:size(amend, 1))
            iterold = changetrackiter(changetrackiter(:, 7) < amend(jj, 7), :); 
            sigiter = significance(iterold(:, 6), :);
            sigiter(:, (amend(jj, 7)):end) = 0;
            significance(iterold(:, 6), :) = sigiter;
        end

        
            
    end

    citation_track = [zeros(length(rda), 1), indtrack, firmtrack, ceil(cda/month_increment), nan(length(rda), 2)];
    indpanel_unique = unique(indtrack);
    for (ii = 1:length(indpanel_unique))
        iter = significance(indtrack == indpanel_unique(ii), :);
        citation_iter = zeros(size(iter, 1), 1);

        for (jj = 1:size(citation_iter, 1))
            if (~any(isnan(iter(jj, :))))
                iterinner = iter(:, find(iter(jj, :) > 0, 1));
                iterinner(jj) = 0;
                citation_iter = citation_iter + (iterinner >= rand(size(iterinner,1), 1));
            end
        end
        citation_track(citation_track(:, 2) == indpanel_unique(ii), 1) = citation_iter;
    end
    citation_track(any(isnan(significance), 2), 1) = nan;    
    citation_track = citation_track(ismember(citation_track(:, 4), 1:(size(innovation_mat, 2)/month_increment - fcitcount + 1 - 5)), :);

    citation_track2 = array2table([(1:size(citation_track, 1))', citation_track(:, 1:4)], ...
                                  'VariableNames', {'wku', 'fcit020', 'CPC', 'Firm', 'Time'});
    citation_track2.count = ones(size(citation_track2, 1), 1); 

    %
    
    citation_track2_orig = p90_rebucket_dev(citation_track2, ...
                                            unique(citation_track2(:, ["Time", "CPC", "Firm"]), 'rows', 'stable'), ...
                                            prctile_inp, method_prc);


    citation_track2 = citation_track2_orig; 
    
    m.citation_track2 = citation_track2; 
   

end