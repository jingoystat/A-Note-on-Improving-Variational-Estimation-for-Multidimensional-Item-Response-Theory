function [new_a, new_b, new_Sigma_theta] = importance_gradient_descent(u, ra, rb, rsigma, mu_i, sig_i, indic, lr, S, M, maxiter, beta_1, beta_2)
% @lr: learning rate
% @S
% @M
% @maxiter
% @eps

    iter = 1;

    new_a = ra;
    new_b = rb;
    new_Sigma_theta = rsigma;
    
    
    [person, item] = size(u);
    [domain, ~] = size(ra);
    
    eps = 0.0001;
    
    v_a = zeros(size(ra));
    v_b = zeros(size(rb));
    v_sigma = zeros(size(rsigma));
    s_a = zeros(size(ra));
    s_b = zeros(size(rb));
    s_sigma = zeros(size(rsigma));

    while iter <= maxiter
        old_a = new_a;
        old_b = new_b;
        old_Sigma_theta = new_Sigma_theta;
        
        %%%%% sampling %%%%%%
        [theta_IS] = importance_sampling(mu_i, sig_i, person, domain, S, M);
        [~, w_tilde] = importance_weights(theta_IS, mu_i, sig_i, old_a, old_b, old_Sigma_theta, u, person, item, S, M);      
        
        [grad_a, grad_b, grad_Sigma_theta] = importance_gradient(u, old_a, old_b, old_Sigma_theta, theta_IS, w_tilde, S, M);        
        
        v_a = beta_1 * v_a + (1 - beta_1) * grad_a;
        v_a = v_a / (1 - beta_1^iter);
        v_b = beta_1 * v_b + (1 - beta_1) * grad_b';
        v_b = v_b / (1 - beta_1^iter);
        v_sigma = beta_1 * v_sigma + (1 - beta_1) * grad_Sigma_theta;
        v_sigma = v_sigma / (1 - beta_1^iter);
        
        s_a = beta_2 * s_a + (1 - beta_2) * grad_a.^2;
        s_a = s_a / (1 - beta_2^iter);
        s_b = beta_2 * s_b + (1 - beta_2) * grad_b'.^2;
        s_b = s_b / (1 - beta_2^iter);
        s_sigma = beta_2 * s_sigma + (1 - beta_2) * grad_Sigma_theta.^2;
        s_sigma = s_sigma / (1 - beta_2^iter);
        
        grad_a_adam = v_a ./ (sqrt(s_a) + 0.001);
        grad_b_adam = v_b ./ (sqrt(s_b) + 0.001);
        grad_sigma_adam = v_sigma ./ (sqrt(s_sigma) + 0.001);
            
        new_a = old_a + lr * grad_a_adam;
        new_a = new_a .* indic;
        new_b = old_b + lr * grad_b_adam;
        new_Sigma_theta = old_Sigma_theta + 0.1 * lr * grad_sigma_adam;
     
        d_temp = sqrt(diag(diag(new_Sigma_theta)));  % standardization matrix
        new_Sigma_theta = inv(d_temp) * new_Sigma_theta * inv(d_temp); %(d_temp\Sigma)/d_temp; %      
        new_Sigma_theta = (new_Sigma_theta + new_Sigma_theta') / 2;
        
        e = max([norm(new_a - old_a), norm(new_b - old_b), norm(new_Sigma_theta - old_Sigma_theta)]);
        if e < eps
            break
        end
   
        iter = iter + 1;
        
    end

       