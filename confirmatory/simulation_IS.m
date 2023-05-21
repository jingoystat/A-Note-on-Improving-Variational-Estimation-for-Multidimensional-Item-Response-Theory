function simulation_IS(person, item, domain, r, within, maxiter, Nrep, S, M)
 
    % learning rates
    lr_list = [0.05]; %[0.1, 0.5, 1];
    N_lr = length(lr_list);

    % a_true, b_true, u, theta_true, Sigma_true
    a_true_list = cell(Nrep,1);
    b_true_list = cell(Nrep,1);
    theta_true_list = cell(Nrep,1);
    Sigma_true_list = cell(Nrep,1);
    data_list = cell(Nrep,1);
    
    % [ra, rb, rsigma, reps, reta, mu_i, sig_i, n]
    ra_list = cell(Nrep,1);
    rb_list = cell(Nrep,1);
    rsigma_list = cell(Nrep,1);
    reta_list = cell(Nrep,1);
    
    
    % [new_a, new_b, new_Sigma, new_mu_i, new_sig_i]
    new_a_list = cell(Nrep,1);
    new_b_list = cell(Nrep,1);
    new_Sigma_list = cell(Nrep,1);
   
    beta_1 = 0.9;
    beta_2 = 0.99;
    
    %%%%%%%%%%%%%% BTW or WITHIN item
    %within = 1;
    %L1_indic = [eye(domain) ones(domain,item-domain)];
    if within == 0
        %%%%%% BTW item model
        % do not consider iteractions
        J=item/domain - 1;
        L1_indic_true = zeros(domain, item-domain);
        for d=1:domain
            L1_indic_true(d,(d-1)*J+1:d*J)=ones(1,J);   % domain*item
        end
        L1_indic_true = [eye(domain) L1_indic_true];
    else   
        %%%%%% WITHIN item model
        % consider interactions'
        if domain == 3
            L1_indic_true = zeros(domain, item-domain);
            L1_indic_true(1,1:8)=1; L1_indic_true(2,9:16)=1; L1_indic_true(3,17:24)=1;
            L1_indic_true(1:2,25:28)=1; L1_indic_true(1,29:32)=1; L1_indic_true(3,29:32)=1;
            L1_indic_true(2:3,33:36)=1; L1_indic_true(:,37:42)=1; 
            L1_indic_true = [eye(domain) L1_indic_true];
        else
            full = [1 0; 0 1; 1 1]';
            L1_indic_true = [eye(domain) full full full full full full];
        end
    end

    for i = 1:Nrep
     
        
        sigL = r;
        sigR = r;
        [a_true, b_true, u, theta_true, Sigma_true] = simulation_2PL(person, item, domain, L1_indic_true, sigL, sigR, within);
        
        %%%%%%%%%%%%% initialization
        [a0, b0, eta0, eps0, Sigma0] = init(u, domain, L1_indic_true);

        %%%%%%%%%%%%% run GVEM function
        [ra, rb, rsigma, reps, reta, mu_i, sig_i, n] = vem_CFA(u, domain, a0, b0, eta0, Sigma0, L1_indic_true);
        
        % a_true, b_true, u, theta_true, Sigma_true
        a_true_list{i} = a_true;
        b_true_list{i} = b_true;
        theta_true_list{i} = theta_true;
        Sigma_true_list{i} = Sigma_true;
        data_list{i} = u;

        % [ra, rb, rsigma, reps, reta, mu_i, sig_i, n]
        ra_list{i} = ra;
        rb_list{i} = rb;
        rsigma_list{i} = rsigma;
        reta_list{i} = reta;

        new_a_list_sub = cell(N_lr, 1);
        new_b_list_sub = cell(N_lr, 1);
        new_sigma_list_sub = cell(N_lr, 1);
        
        for j = 1:N_lr
            
            lr = lr_list(j);
            
            %%%% importance sampling %%%%
            [new_a, new_b, new_sigma, corr_list, a_bias_list, b_bias_list, corr_bias_list, lb_list] = importance_gradient_descent(u, a_true, b_true, Sigma_true, ra, rb, rsigma, mu_i, sig_i, L1_indic_true, lr, S, M, maxiter, beta_1, beta_2);
            % [new_a, new_b, new_Sigma]
            new_a_list_sub{j} = new_a;
            new_b_list_sub{j} = new_b;
            new_sigma_list_sub{j} = new_sigma;
        end
        
        new_a_list{i} = new_a_list_sub;
        new_b_list{i} = new_b_list_sub;
        new_Sigma_list{i} = new_sigma_list_sub;
    end 
    
    if within == 0
        model = 'between';
    else
        model = 'within';
    end 
    
    save("GVEM_IS_"+person+"_"+item+"_"+domain+"_"+r+"_"+model+"_local.mat");
