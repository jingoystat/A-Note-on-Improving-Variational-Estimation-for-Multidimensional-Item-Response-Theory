function [theta_IS] = importance_sampling(mu_i, sig_i, person, domain, S, M)
    
    theta_IS = arrayfun(@(i) reshape(mvnrnd(mu_i(:, i), squeeze(sig_i(:,:,i)), S*M), [S, M, 1, domain]), [1:person], 'UniformOutput',false);
    theta_IS = cat(3, theta_IS{:});
end

