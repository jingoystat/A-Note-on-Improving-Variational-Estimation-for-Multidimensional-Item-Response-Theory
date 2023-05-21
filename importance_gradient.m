function [grad_a, grad_b, grad_Sigma_theta] = importance_gradient(y, a, b, Sigma_theta, theta_IS, w_tilde, S, M)
    [person, item] = size(y);
    [domain, ~] = size(a);
    
    tmp = squeeze(sum(bsxfun(@times, reshape(theta_IS, [S, M, person, domain, 1]), reshape(a, [1, 1, 1, domain, item])), 4)); %%%% S, M, person, item
    partial = reshape(y, [person, 1, 1, item]) - 1 + 1 / (1 +  exp(permute(tmp, [3,1,2,4]) - reshape(b, [1,1,1,item]))); % person, S, M, item
    partial = permute(partial, [2,3,1,4]); % S, M, person, item
    grad_a = bsxfun(@times, reshape(partial, [S, M, person, 1, item]), reshape(theta_IS, [S, M, person, domain, 1]));
    w_tilde = permute(w_tilde, [2,3,1]); % S * M * person
    
    grad_a = squeeze(sum(bsxfun(@times, reshape(w_tilde, [S, M, person, 1, 1]), grad_a), 2));
    grad_a = squeeze(mean(grad_a, [1,2]));
    
    grad_b = squeeze(sum(bsxfun(@times, reshape(w_tilde, [S, M, person, 1]), partial), 2));
    grad_b = -squeeze(mean(grad_b, [1,2]));
    
    tmp_sigma = bsxfun(@times, reshape(theta_IS, [S, M, person, domain, 1]), reshape(theta_IS, [S, M, person, 1, domain])); % S, M, person, domain, domain
    grad_sigma = bsxfun(@times, reshape(w_tilde, [S, M, person, 1, 1]),(reshape(Sigma_theta, [1,1,1,domain, domain]) - tmp_sigma)/2);
    grad_sigma = squeeze(sum(grad_sigma, 2));
    grad_Sigma_theta = squeeze(mean(grad_sigma, [1,2]));
    

    
    
