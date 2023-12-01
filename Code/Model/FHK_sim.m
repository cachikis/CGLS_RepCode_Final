%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: FHK_sim.m
% Author: Craig A. Chikis
% Date: 08/10/2021
% Note(s):
% MA-MFA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function m = FHK_sim(m)


    innovation_mat = m.innovation_mat;
    scale_time = m.scale_time;
    leader_innovL = m.leader_innovL;
    laggard_innovL = m.laggard_innovL;
    entrant_innovL = m.entrant_innovL;
    leader_innovF = m.leader_innovF;
    laggard_innovF = m.laggard_innovF;
    entrant_innovF = m.entrant_innovF;
    lambda = m.lambda; 
    ageL = m.ageL;
    ageF = m.ageF;
    nu_s = m.nu_s; 
    kappa = m.kappa; 
    s_max = m.s_max; 
    LMS_bar_s = m.LMS_bar_s; 



	width_param = round(5*12/scale_time);

	stateL = innovation_mat(:, :, 1, 2);
	stateF = innovation_mat(:, :, 4, 2); 
	state_all = [stateL; stateF];



	age_all = [ageL; ageF];
	age_all = age_all(:, (end-(width_param-1)):end);

	leader_innov = [leader_innovL; leader_innovF];
	laggard_innov = [laggard_innovL; laggard_innovF];
	entrant_innov = [entrant_innovL; entrant_innovF];

	leader_innov = leader_innov(:, (end-(width_param-1)):end);
	laggard_innov = laggard_innov(:, (end-(width_param-1)):end);
	entrant_innov = entrant_innov(:, (end-(width_param-1)):end);

	state_all = state_all(:, (end-(width_param-1)):end);

	C_firms = zeros(size(age_all,1), 1);
	N_firms = zeros(size(age_all,1), 1);

	for (ii = 1:size(C_firms,1))
		C_firms(ii) = ~ismember(0, age_all(ii, :));
		N_firms(ii) = ismember(0, age_all(ii,:));
	end

	C_firms = find(C_firms == 1);
	N_firms = find(N_firms == 1);



	C_panel = state_all(C_firms, :);
	N_panel = state_all(N_firms, :);

	N_age = age_all(N_firms, :);

	N_where = zeros(size(N_age,1), 1);
	for (ii = 1:size(N_age,1))
		N_where(ii) = find(N_age(ii,:) == 0,1);
	end

	leader_innov_C = leader_innov(C_firms, :);
	laggard_innov_C = laggard_innov(C_firms, :);
	entrant_innov_C = entrant_innov(C_firms, :);

	leader_innov_N = leader_innov(N_firms, :);
	laggard_innov_N = laggard_innov(N_firms, :);
	entrant_innov_N = entrant_innov(N_firms, :);

	if (kappa >= 9999)
		weightL = ones(1, s_max);
		weightF = zeros(1, s_max);
		weight0 = .5;

		weights = [flip(weightF), weight0, weightL];
		act = weights;
	else
		weights = [flip(1./(1+nu_s.^(-kappa))), 1./(1+nu_s(2:end).^(kappa))];
		act = weights;
	end



	theta_tm1 = weights(s_max+1+state_all(:,1))';
			

	theta_t = weights(s_max+1+state_all(:,end))';

	Dtheta_t = theta_t - theta_tm1;


	steps_C = zeros(size(C_panel));
	steps_N = zeros(size(N_panel));
	steps_all_ind = zeros(size([steps_C; steps_N]));
	steps_all_ind_alt = zeros(size(steps_all_ind));
	steps_all_ind_alt(:,1) = lambda.^state_all(:,1);

	steps_N(:, 1) = 1;

	stateL_all = state_all(1:size(ageL,1), :);
	stateF_all = state_all((size(ageL,1)+1):end, :);

	ageL_all = age_all(1:size(ageL,1), :);
	ageF_all = age_all((size(ageL,1)+1):end, :);

	stateL_all = reshape(stateL_all', 1, size(stateL_all,2), size(stateL_all,1));
	stateF_all = reshape(stateF_all', 1, size(stateF_all,2), size(stateF_all,1));

	ageL_all = reshape(ageL_all', 1, size(ageL_all,2), size(ageL_all,1));
	ageF_all = reshape(ageF_all', 1, size(ageF_all,2), size(ageF_all,1));

	state_all_new = zeros(2, size(stateL_all,2), size(stateL_all,3));
	state_all_new(1,:, :) = stateL_all(1, :, :);
	state_all_new(2,:, :) = stateF_all(1, :, :);

	age_all_new = zeros(2, size(ageL_all,2), size(ageL_all,3));
	age_all_new(1,:, :) = ageL_all(1, :, :);
	age_all_new(2,:, :) = ageF_all(1, :, :);

	prod_track = nan(size(state_all_new(:,:,1)));
	prod_track_outer = nan(size(state_all_new));
	prod_track_outer_entering = nan(size(state_all_new));

	for (ii = 1:size(state_all_new, 3))
		s_iter = state_all_new(:, :, ii);
		age_iter = age_all_new(:, :, ii);

		find_laggard = find(s_iter(:,1) <= 0);
		find_leader = find(s_iter(:,1) >= 0);

		find_entrant1 = ismember(0, age_iter(1,2:end));
		find_entrant2 = ismember(0, age_iter(2,2:end));

        
        if (abs(s_iter(1, 1)) > LMS_bar_s)
            prod_track(find_laggard, 1) = abs(s_iter(find_laggard, 1))+0;
            prod_track(find_leader, 1) = abs(s_iter(find_laggard, 1))+LMS_bar_s;
        else
            prod_track(find_laggard, 1) = abs(s_iter(find_laggard,1))+0;
            prod_track(find_leader, 1) = (abs(s_iter(find_laggard,1))+0)*2; 
        end
	

		for (jj = 2:size(state_all_new,2))
            if (abs(s_iter(1, jj)) > LMS_bar_s)
				change = max(s_iter(:, jj) - s_iter(:, jj-1), 0).*(s_iter(:, jj) > 0);
                prod_track(:, jj) = prod_track(:, jj-1) + max(change, [], 'all');
            else	
				change = max(s_iter(:, jj) - s_iter(:, jj-1), 0);
                prod_track(:, jj) = prod_track(:, jj-1) + change;
            end


		end

		if (find_entrant1)
			prod_track_outer_entering(1,:,ii) = prod_track(1,:);
		else
			prod_track_outer(1,:,ii) = prod_track(1,:);
		end

		if (find_entrant2)
			prod_track_outer_entering(2,:,ii) =prod_track(2,:);
		else
			prod_track_outer(2,:,ii) = prod_track(2,:);
		end

	end
	steps = prod_track_outer;
	steps_enter = prod_track_outer_entering;

	prod_line = lambda.^prod_track_outer;
	prod_line_N = lambda.^prod_track_outer_entering;

	prod_track_outer(1,:,:) = weights(s_max+1+state_all_new(1,:,:)).*prod_line(1,:,:);
	prod_track_outer(2,:,:) = weights(s_max+1+state_all_new(2,:,:)).*prod_line(2,:,:);

	prod_track_outer_entering(1,:,:) = weights(s_max+1+state_all_new(1,:,:)).*prod_line_N(1,:,:);
	prod_track_outer_entering(2,:,:) = weights(s_max+1+state_all_new(2,:,:)).*prod_line_N(2,:,:);


	line1 = reshape(nansum(nansum([prod_track_outer(:,1,:), prod_track_outer_entering(:,1,:)], 1),2), ...
					size(ageL,1),1);
	line2 = reshape(nansum(nansum([prod_track_outer(:,end,:), prod_track_outer_entering(:,end,:)], 1),2), ...
					size(ageL,1),1);

	prod_growth_mean = mean(log(max(prod_track_outer(:, end, :), [], 1)) - log(max(prod_track_outer(:, 1, :), [], 1)));



	WITHIN_vec = reshape(nansum(act(s_max+1+state_all_new(:,1,:)).*((prod_line(:,end,:)) - (prod_line(:,1,:))), 1), ...
						 size(ageL,1),1);
	
	BETWEEN_vec = reshape(nansum(((prod_line(:,1,:)) - (reshape(line1, 1, 1, length(line1)))).*... 
						 	  (act(s_max+1+state_all_new(:,end,:)) - act(s_max+1+state_all_new(:,1,:))),1), ...
						  size(ageL,1), 1);
	CROSS_vec = reshape(nansum(((prod_line(:,end,:)) - (prod_line(:,1,:))).*...
							(act(s_max+1+state_all_new(:,end,:)) - act(s_max+1+state_all_new(:,1,:))),1), ...
						size(ageL,1),1);

	ENTER_vec = reshape(nansum(act(s_max+1+state_all_new(:,end,:)).*...
						((prod_line_N(:,end,:)) - (reshape(line1, 1, 1, length(line1)))),1), ...
						size(ageL,1),1);

	EXIT_vec = reshape(nansum(act(s_max+1+state_all_new(:,1,:)).*...
						((prod_line_N(:,1,:)) - (reshape(line1, 1, 1, length(line1)))),1), ...
						size(ageL,1),1);

	check = [((line2) - (line1)), ...
			 -WITHIN_vec, -BETWEEN_vec, -CROSS_vec, -ENTER_vec, EXIT_vec];

	WITHIN_vec_log = reshape(nansum(act(s_max+1+state_all_new(:,1,:)).*(log(prod_line(:,end,:)) - log(prod_line(:,1,:))), 1), ...
						 size(ageL,1),1);
	
	BETWEEN_vec_log = reshape(nansum((log(prod_line(:,1,:)) - log(reshape(line1, 1, 1, length(line1)))).*... 
						 	  (act(s_max+1+state_all_new(:,end,:)) - act(s_max+1+state_all_new(:,1,:))),1), ...
						  size(ageL,1), 1);
	CROSS_vec_log = reshape(nansum((log(prod_line(:,end,:)) - log(prod_line(:,1,:))).*...
							(act(s_max+1+state_all_new(:,end,:)) - act(s_max+1+state_all_new(:,1,:))),1), ...
						size(ageL,1),1);

	ENTER_vec_log = reshape(nansum(act(s_max+1+state_all_new(:,end,:)).*...
						(log(prod_line_N(:,end,:)) - log(reshape(line1, 1, 1, length(line1)))),1), ...
						size(ageL,1),1);

	EXIT_vec_log = reshape(nansum(act(s_max+1+state_all_new(:,1,:)).*...
						(log(prod_line_N(:,1,:)) - log(reshape(line1, 1, 1, length(line1)))),1), ...
						size(ageL,1),1);

	check_log = [(log(line2) - log(line1)), ...
			 -WITHIN_vec_log, -BETWEEN_vec_log, -CROSS_vec_log, -ENTER_vec_log, EXIT_vec_log];



	[check_sort, check_idx] = sort(abs(sum(check,2)), 'descend');
	tmp	= find(abs(sum(check,2)) > 1e-6);

	[check_sort_log, check_idx_log] = sort(abs(sum(check_log,2)), 'descend');
	tmp_log	= find(abs(sum(check_sort_log,2)) > 1e-6);

	if (~(kappa >= 9999))
		WITHIN = WITHIN_vec;
		BETWEEN = BETWEEN_vec;
		CROSS = CROSS_vec;
		ENTRY = ENTER_vec;
		EXIT = EXIT_vec;
	else
		WITHIN = WITHIN_vec_log;
		BETWEEN = BETWEEN_vec_log;
		CROSS = CROSS_vec_log;
		ENTRY = ENTER_vec_log;
		EXIT = EXIT_vec_log;
	end

	% Equal weighted mean across industries (log preferences ==> equal revenue weights)
	WITHIN_v = WITHIN ./ (WITHIN+BETWEEN+CROSS+ENTRY-EXIT);
	WITHIN_v = WITHIN_v(~isinf(WITHIN_v)); 

	BETWEEN_v = BETWEEN ./ (WITHIN+BETWEEN+CROSS+ENTRY-EXIT);
	BETWEEN_v = BETWEEN_v(~isinf(BETWEEN_v)); 

	CROSS_v = CROSS ./ (WITHIN+BETWEEN+CROSS+ENTRY-EXIT);
	CROSS_v = CROSS_v(~isinf(CROSS_v)); 

	ENTRY_v = ENTRY ./ (WITHIN+BETWEEN+CROSS+ENTRY-EXIT);
	ENTRY_v = ENTRY_v(~isinf(ENTRY_v)); 

	EXIT_v = EXIT ./ (WITHIN+BETWEEN+CROSS+ENTRY-EXIT);
	EXIT_v = EXIT_v(~isinf(EXIT_v)); 

    m.WITHIN_out = sum(WITHIN)/sum(WITHIN+BETWEEN+CROSS+ENTRY-EXIT); 
    m.BETWEEN_out = sum(BETWEEN)/sum(WITHIN+BETWEEN+CROSS+ENTRY-EXIT); 
    m.CROSS_out = sum(CROSS)/sum(WITHIN+BETWEEN+CROSS+ENTRY-EXIT); 
    m.ENTRY_out = sum(ENTRY)/sum(WITHIN+BETWEEN+CROSS+ENTRY-EXIT); 
    m.EXIT_out = sum(EXIT)/sum(WITHIN+BETWEEN+CROSS+ENTRY-EXIT); 

end
		


