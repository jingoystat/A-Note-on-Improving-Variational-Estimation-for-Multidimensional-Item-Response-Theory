function simulation_IS_parallel(person, item, domain, r, within, maxiter, Nrep, S, M)

    lr_list = [0.5, 0.1, 0.05];
    N_lr = length(lr_list);
    a_true_list = cell(Nrep,1);
    b_true_list = cell(Nrep,1);
    theta_true_list = cell(Nrep,1);
    Sigma_true_list = cell(Nrep,1);
    data_list = cell(Nrep,1);
    
    ra_list = cell(Nrep,1);
    rb_list = cell(Nrep,1);
    rsigma_list = cell(Nrep,1);
    reta_list = cell(Nrep,1);
    
    
    new_a_list = cell(Nrep,1);
    new_b_list = cell(Nrep,1);
    new_Sigma_list = cell(Nrep,1);
   
    beta_1 = 0.9;
    beta_2 = 0.999;
    
    converge_list = zeros(Nrep, 1);
    success_list = zeros(Nrep, 1);
 
    
    time_list = zeros(Nrep, 2);
    
    if within == 0
        J=item/domain - 1;
        indic = zeros(domain, item-domain);
        for d=1:domain
            indic(d,(d-1)*J+1:d*J)=ones(1,J);  
        end
        indic = [eye(domain) indic];
    else   
       
        if domain == 3
            indic = zeros(domain, item-domain);
            indic(1,1:5)=1; indic(2,6:10)=1; indic(3,11:15)=1;
            indic(1:2,16:21)=1; 
            indic(1,22:27)=1; indic(3,22:27)=1;
            indic(2:3,28:33)=1; 
            indic(:,34:45)=1;
        elseif domain == 2
            full = [1 0; 0 1; 1 1]';
            indic = [];
            for j = 1:item/3
                indic = [indic full];
            end
        else
            indic = zeros(domain, item-domain);
            indic(1, 1:2) = 1; indic(2, 3:4) = 1; indic(3, 5:6) = 1; indic(4, 7:8) = 1; indic(5, 9:10) = 1;
            indic(1:2, 11:12) = 1; indic([1,3], 13:14) = 1; indic([1,4], 15:16) = 1; indic([1,5], 17:18) = 1; 
            indic(2:3, 19:20) = 1; indic([2,4], 21:22) = 1; indic([2,5], 23:24) = 1; 
            indic(3:4, 25:26) = 1; indic([3,5], 27:28) = 1; 
            indic(4:5, 29:30) = 1;
            indic([1,2,3], 31:32) = 1; indic([1,2,4], 33:34) = 1; indic([1,2,5], 35:36) = 1; 
            indic([1,3,4], 37:38) = 1; indic([1,3,5], 39:40) = 1; indic([1,4,5], 41:42) = 1;
            indic([2,3,4], 43:44) = 1; indic([2,3,5], 45:46) = 1; indic([2,4,5], 47:48) = 1;
            indic([3,4,5], 49:50) = 1;
            indic = [eye(domain) indic]; 
        end
    end

    if within == 0
        model = 'between';
    else
        model = 'within';
    end 

  
    parfor i = 1:Nrep
      
        if r == "high"
            sigL = 0.5;
            sigR = 0.7;
        else
            sigL = 0.1;
            sigR = 0.3;
        end

        [a_true, b_true, u, theta_true, Sigma_true] = simulation_2PL(person, item, domain, indic, sigL, sigR, within);
        
        writematrix(u, "data_"+person+"_"+item+"_"+domain+"_"+r+"_"+model+"_rep"+i+".csv");

        tStart = cputime;
        %%%%%%%%%%%%% initialization
        [a0, b0, eta0, eps0, ~] = init(u, domain, ones(domain,item));

        try
        %%%%%%%%%%%%% run GVEM function
        [ra, rb, rsigma, reps, reta, mu_i, sig_i, n, converge]=vem_EFAro(u, domain, a0, b0, eta0);
        
        tGVEM = cputime - tStart;
      
        a_true_list{i} = a_true;
        b_true_list{i} = b_true;
        theta_true_list{i} = theta_true;
        Sigma_true_list{i} = Sigma_true;
        data_list{i} = u;
        converge_list(i) = converge;

        ra_list{i} = ra;
        rb_list{i} = rb;
        rsigma_list{i} = rsigma;
        reta_list{i} = reta;
         
        best_lb = -Inf;
        tStart = cputime;  
        
        theta_eval = importance_sampling(mu_i, sig_i, person, domain, S, M);
        best_lr = 0;
        for j = 1:N_lr
            lr = lr_list(j); 
            [new_a, new_b, new_sigma] = importance_gradient_descent(u, ra, rb, rsigma, mu_i, sig_i, lr, S, M, maxiter, beta_1, beta_2);
           
            [w, ~] = importance_weights(theta_eval, mu_i, sig_i, new_a, new_b, new_sigma, u, person, item, S, M);
            lb_val = lower_bound(w);
            
            if lb_val > best_lb
                new_a_list{i} = new_a;
                new_b_list{i} = new_b;
                new_Sigma_list{i} = new_sigma;
                best_lb = lb_val;
                best_lr = lr;
            end
        end
       
        tIS = cputime - tStart;
        time_list(i,:) = [tGVEM, tIS];
        success_list(i) = 1;
        catch
            success_list(i) = 0;
        end
    end
    
    
   
    
    save("GVEM_IS_"+person+"_"+item+"_"+domain+"_"+r+"_"+model+".mat");
