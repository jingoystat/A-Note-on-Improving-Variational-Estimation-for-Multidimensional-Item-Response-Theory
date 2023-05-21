function [new_a, new_b, new_Sigma_theta] = importance_stochastic_gradient_descent(u, ra, rb, rsigma, mu_i, sig_i, lr, S, M, maxiter)
% @lr: learning rate
% @S
% @M
% @maxiter
% @eps

    iter = 0;


    new_a = ra;
    new_b = rb;
    new_Sigma_theta = rsigma;
    
    [person, ~] = size(u);
    [domain, ~] = size(ra);
    
    eps = 0.001;
    n_stochastic = 200;

    while iter < maxiter
        old_a_before = new_a;
        old_b_before = new_b;
        old_Sigma_theta_before = new_Sigma_theta;
        
        %%%%% sampling %%%%%%
        theta_IS = importance_sampling(mu_i, sig_i, person, domain, S, M);
        
        %%%%% gradient descent %%%%%
        iter_grad = 1;
        while iter_grad < 20
            prop = n_stochastic / person;
            index = binornd(1, prop, person);
            y_sto = 
            mu_i_sto
            sig_i_sto = 
            
            old_a = new_a;
            old_b = new_b;
            old_Sigma_theta = new_Sigma_theta;
        
            [grad_a, grad_b, grad_Sigma_theta] = importance_gradient(u, old_a, old_b, old_Sigma_theta, mu_i, sig_i, theta_IS, S, M);
            new_a = old_a + lr * grad_a';
            new_b = old_b + lr * grad_b';
            new_Sigma_theta = old_Sigma_theta + lr * grad_Sigma_theta;
            d_temp = sqrt(diag(diag(new_Sigma_theta)));  % standardization matrix
            new_Sigma_theta = inv(d_temp) * new_Sigma_theta * inv(d_temp); %(d_temp\Sigma)/d_temp; %      
            new_Sigma_theta = (new_Sigma_theta + new_Sigma_theta') / 2;
            
            % new_Sigma_theta
            [~,flag] = chol(new_Sigma_theta);
            if flag ~= 0
                new_Sigma_theta = old_Sigma_theta;
            end
            %flag, new_Sigma_theta
            
            iter_grad = iter_grad + 1;
            e_grad = max([norm(grad_a), norm(grad_b), norm(grad_Sigma_theta)]);
            if e_grad < eps
                break
            end
        end
       
        iter = iter + 1;
        e = max([norm(old_a_before - new_a), norm(old_b_before - new_b), norm(old_Sigma_theta_before - new_Sigma_theta)]);
        if e < eps
            break
        end
    end
