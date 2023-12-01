%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Filename: bucket_quantile.m
% Author: Craig A. Chikis
% Date: 08/27/2022
% Note(s):
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function table_in = bucket_quantile(table_in, col_t, nbucket)

  
    col_tm1 = strcat(col_t, "m1");
    col_t1 = strcat(col_t, "1"); 

    uniqueyears = unique(table_in.Time); 
    table_in.bucket = nan(size(table_in, 1), 1);
    table_in.bucket_tm1 = nan(size(table_in, 1), 1);
    table_in.bucket_t1 = nan(size(table_in, 1), 1);
    for (ii = 1:length(uniqueyears))
        iter_t = table_in(table_in.Time == uniqueyears(ii) & ...
                          ~isnan(table2array(table_in(:, col_t))), [col_t, "bucket"]);
        iter_tm1 = table_in(table_in.Time == uniqueyears(ii) & ...
                            ~isnan(table2array(table_in(:, col_tm1))), [col_tm1, "bucket_tm1"]);
        iter_t1 = table_in(table_in.Time == uniqueyears(ii) & ...
                           ~isnan(table2array(table_in(:, col_t1))), [col_t1, "bucket_t1"]);


        [~, idx_t] = sortrows(iter_t, col_t, 'descend');
        [~, idx_tm1] = sortrows(iter_tm1, col_tm1, 'descend');
        [~, idx_t1] = sortrows(iter_t1, col_t1, 'descend'); 

        ct_t = floor(length(idx_t)/nbucket);
        ct_tm1 = floor(length(idx_tm1)/nbucket);
        ct_t1 = floor(length(idx_t1)/nbucket); 
        st = 0; 
        for (bucketiter = 1:nbucket)
            if (bucketiter < nbucket)
                iter_t.bucket(idx_t((st+1):(nbucket*ct_t))) = bucketiter;
            else
                iter_t.bucket(idx_t((st+1):end)) = bucketiter;
            end


            st = st + ct_t; 
        end

        iter_t = iter_t(:, "bucket"); 

   
        st = 0; 
        for (bucketiter = 1:nbucket)
            if (bucketiter < nbucket)
                iter_tm1.bucket_tm1(idx_tm1((st+1):(nbucket*ct_tm1))) = bucketiter;
            else
                iter_tm1.bucket_tm1(idx_tm1((st+1):end)) = bucketiter;
            end
            st = st + ct_tm1; 
        end
        iter_tm1 = iter_tm1(:, "bucket_tm1"); 

        st = 0; 
        for (bucketiter = 1:nbucket)
            if (bucketiter < nbucket)
                iter_t1.bucket_t1(idx_t1((st+1):(nbucket*ct_t1))) = bucketiter;
            else
                iter_t1.bucket_t1(idx_t1((st+1):end)) = bucketiter;
            end
            iter_t1.bucket_t1(idx_t1((st+1):(nbucket*ct_t1))) = bucketiter;
            st = st + ct_t1; 
        end
        iter_t1 = iter_t1(:, "bucket_t1"); 

        table_in.bucket(table_in.Time == uniqueyears(ii) & ...
                        ~isnan(table2array(table_in(:, col_t)))) = iter_t.bucket;
                    
        table_in.bucket_tm1(table_in.Time == uniqueyears(ii) & ...
                            ~isnan(table2array(table_in(:, col_tm1)))) = iter_tm1.bucket_tm1;
        
        table_in.bucket_t1(table_in.Time == uniqueyears(ii) & ...
                           ~isnan(table2array(table_in(:, col_t1)))) = iter_t1.bucket_t1;

    end



end