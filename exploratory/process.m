
% a_bias = zeros(Nrep, domain);
a_bias = zeros(Nrep, 1);
b_bias = zeros(Nrep, 1);
corr_bias = zeros(Nrep, 1);

% ra_bias = zeros(Nrep, domain);
ra_bias = zeros(Nrep, 1);
rb_bias = zeros(Nrep, 1);
rcorr_bias = zeros(Nrep, 1);


% a_rmse = zeros(Nrep, domain);
a_rmse = zeros(Nrep, 1);
b_rmse = zeros(Nrep, 1);
corr_rmse = zeros(Nrep, 1);

% ra_rmse = zeros(Nrep, domain);
ra_rmse = zeros(Nrep, 1);
rb_rmse = zeros(Nrep, 1);
rcorr_rmse = zeros(Nrep, 1);


N_converge = 0;




for i = 1:Nrep
    
    if converge_list(i) == 0  && success_list(i) == 1 && isreal(ra_list{i}) && isreal(new_a_list{i})
            %i
            N_converge = N_converge + 1;
            [ra, rsigma] = rotate_A_Sig(ra_list{i}, rsigma_list{i}, a_true_list{i});
            % ra_bias(i,:) = sum(ra_list{i} - a_true_list{i}, 2) ./ sum(L1_indic_true == 1, 2);
            ra_bias(i) = mean(ra - a_true_list{i}, 'all');%  ./ sum(L1_indic_true == 1, 'all');
            rb_bias(i) = mean(rb_list{i}- b_true_list{i}, 'all');

            [new_a, new_sigma] = rotate_A_Sig(new_a_list{i}, new_Sigma_list{i}, a_true_list{i});
            % a_bias(i,:) = sum(new_a_list{i} - a_true_list{i}, 2) ./ sum(L1_indic_true == 1, 2);
            a_bias(i) = mean(new_a - a_true_list{i}, 'all');% ./ sum(L1_indic_true == 1, 'all');
            b_bias(i) = mean(new_b_list{i}- b_true_list{i}, 'all');

            if domain == 2
                rcorr_bias(i) = rsigma_list{i}(1, 2) - Sigma_true_list{i}(1, 2);     
                corr_bias(i) = new_sigma(1, 2) - Sigma_true_list{i}(1, 2);  
            else
                rcorr_bias(i) = rsigma(1, 2) - Sigma_true_list{i}(1, 2) + rsigma(1, 3) - Sigma_true_list{i}(1, 3); % + rsigma(2, 3) - Sigma_true_list{i}(2, 3);
                rcorr_bias(i) = rcorr_bias(i) / 2;
                corr_bias(i) = new_sigma(1, 2) - Sigma_true_list{i}(1, 2) + new_sigma(1, 3) - Sigma_true_list{i}(1, 3); % + new_sigma(2, 3) - Sigma_true_list{i}(2, 3);
                corr_bias(i) = corr_bias(i) / 2;
            end


            % a_rmse(i, :) = sqrt(sum((new_a_list{i} - a_true_list{i}).^2, 2) ./ sum(L1_indic_true == 1, 2)) ;
            a_rmse(i) = sqrt(mean((new_a- a_true_list{i}).^2, 'all'));% ./ sum(L1_indic_true == 1, 'all')) ;
            b_rmse(i) = sqrt(mean((new_b_list{i}- b_true_list{i}).^2, 'all'));
            % ra_rmse(i, :) = sqrt(sum((ra_list{i} - a_true_list{i}).^2, 2)  ./ sum(L1_indic_true == 1, 2));
            ra_rmse(i, :) = sqrt(mean((ra - a_true_list{i}).^2, 'all')); %  ./ sum(L1_indic_true == 1, 'all'));
            rb_rmse(i) = sqrt(mean((rb_list{i}- b_true_list{i}).^2, 'all'));

            if domain == 2
                corr_rmse(i) = abs(new_sigma(1, 2) - Sigma_true_list{i}(1, 2));     
            else
                corr_rmse(i) = (new_sigma(1, 2) - Sigma_true_list{i}(1, 2))^2 + (new_sigma(1, 3) - Sigma_true_list{i}(1, 3))^2; % + (new_sigma(2, 3) - Sigma_true_list{i}(2, 3))^2;
                corr_rmse(i) = sqrt(corr_rmse(i) / 2);
            end


            if domain == 2
                rcorr_rmse(i) = abs(rsigma(1, 2) - Sigma_true_list{i}(1, 2));     
            else
                rcorr_rmse(i) = (rsigma(1, 2) - Sigma_true_list{i}(1, 2))^2 + (rsigma(1, 3) - Sigma_true_list{i}(1, 3))^2; % + (rsigma(2, 3) - Sigma_true_list{i}(2, 3))^2;
                rcorr_rmse(i) = sqrt(rcorr_rmse(i) / 2);
            end
    end
    
end


times = mean(time_list,1);
N_converge

res1 = [mean(ra_bias), mean(a_bias), mean(rb_bias), mean(b_bias), mean(rcorr_bias), mean(corr_bias), times(1), times(2)] * Nrep / N_converge 
res2 = [mean(ra_rmse), mean(a_rmse), mean(rb_rmse), mean(b_rmse), mean(rcorr_rmse), mean(corr_rmse)] * Nrep / N_converge 

