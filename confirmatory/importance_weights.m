function [w, w_tilde] = importance_weights(theta_IS, mu_i, sig_i, a, b, Sigma_theta, y, person, item, S, M)
    

    domain = size(a, 1);
    
    kron_b = reshape(kron(ones(person,1),b), [1,1,person,item]); 
    p_tmp = squeeze(-sum(bsxfun(@times, reshape(theta_IS, [S, M, person, domain, 1]), reshape(a, [1,1,1,domain, item])), 4)) + kron_b;
    p_item = 1 / (1 + exp(p_tmp));
    
    part1 = sum(bsxfun(@times, reshape(y, [1,1,person,item]), log(p_item)), 4);
    part2 = sum(bsxfun(@times, reshape(1-y, [1,1,person,item]), log(1-p_item)), 4);
    part1 = permute(part1, [3,1,2]);
    part2 = permute(part2, [3,1,2]);
    
    log_det_sigma_theta = log(det(Sigma_theta));
    inv_sigma_theta = inv(Sigma_theta);
    log_det_sigma_i_fun = @(i) log(det(sig_i(:,:,i)));
    log_det_sigma_i = arrayfun(log_det_sigma_i_fun, [1:person], 'UniformOutput',false);
    inv_sigma_i_fun = @(i) inv(sig_i(:,:,i));
    inv_sigma_i_cell = arrayfun(inv_sigma_i_fun, [1:person], 'UniformOutput',false);
    inv_sigma_i = zeros(person, domain, domain);
    for i = 1:person
        inv_sigma_i(i,:,:) = inv_sigma_i_cell{i};
    end
    
    part3 = -squeeze(sum(bsxfun(@times, squeeze(sum(bsxfun(@times, reshape(theta_IS,[S, M, person, domain, 1]), reshape(inv_sigma_theta, [1,1,1,domain,domain])), 4)), theta_IS), 4)); % S, M, person
    part3 = part3/2 - log_det_sigma_theta/2;
    part3 = permute(part3, [3,1,2]);
    
   
    mu_i = permute(mu_i, [2,1]); 
    theta_tmp = theta_IS - reshape(mu_i, [1,1,person,domain]);
    q = squeeze(sum(bsxfun(@times, squeeze(sum(bsxfun(@times, reshape(theta_tmp, [S, M, person, domain, 1]), reshape(inv_sigma_i, [1,1,person,domain,domain])), 4)), theta_tmp), 4));

    q = -q/2 - reshape(cell2mat(log_det_sigma_i)/2, [1, 1, person]);
    q = permute(q, [3,1,2]);
    
    
    p = part1 + part2 + part3;
    w = exp(p - q);
    w_tilde = w ./ sum(w, 3);
    
   
    
    