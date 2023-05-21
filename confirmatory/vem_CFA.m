%GVEM algorithm for 2PL model - CFA. 
function [ra,rb,rsigma,reps,reta,mu_i,sig_i,n,converge]=vem_CFA(u,domain,new_a,new_b,eta,Sigma,indic)
    [person,item]=size(u);
    converge = 0;
    xi = eta;  
    MU = zeros(domain, person);
    SIGMA = zeros(domain, domain, person);
    n = 0;
  
    
    par_Sigma = Sigma;
    
    new_a = new_a .* indic;


    while converge == 0 && rank(Sigma) == domain  

        par_Sigma=Sigma; 
        
        Spart=zeros(domain,domain);
        
        for i=1:person 
            sigma_part = zeros(domain,domain);
            mu_part = zeros(domain,1);
            
            for j = 1:item
                sigma_part = sigma_part + eta(i,j) * new_a(:,j) * new_a(:,j)';    
                mu_part = mu_part + (2 * eta(i,j) * new_b(j) + u(i,j) - 0.5) * new_a(:,j);     
            end
            
            sigmahat = inv(inv(Sigma) + 2 * sigma_part);   
            
            muhat= sigmahat * mu_part; 
            
            SIGMA(:,:,i) = sigmahat;
            
            MU(:,i) = muhat;
            xi(i,:) = sqrt(new_b .^ 2 - 2 * new_b .* (new_a' * muhat)' + diag(new_a' * (sigmahat + muhat * muhat') * new_a)');
            Spart = Spart + sigmahat + muhat * muhat';
        end

        eta = (exp(xi) ./ (1 + exp(xi)) - 0.5) ./ (2 * xi);  
        for i = 1:person
            for j = 1:item
                if abs(xi(i,j)) < 0.01
                    eta(i,j) = 0.125;  
                end
            end    
        end

        Sigma = Spart / person;
        d_temp = sqrt(diag(diag(Sigma)));  
        Sigma = inv(d_temp) * Sigma * inv(d_temp);            

       
        par_b=new_b;
        b_part = zeros(person,item);
        for i=1:person
            b_part(i,:)=(new_a'*MU(:,i))';  
        end
        new_b=sum(0.5-u+2*eta.*b_part)./sum(2*eta);  
        clear i              

       
        par_a=new_a;
        for j=1:item
            a_nu=0;a_de=0;
            Ind=indic(:,j);
            iind=find(Ind==1);
            for i=1:person
                sigma=SIGMA(:,:,i);
                sigma=sigma(iind,iind);
                mu=MU(iind,i);
                a_de=a_de+eta(i,j)*sigma+eta(i,j)*(mu*mu');
                a_nu=a_nu+(u(i,j)-0.5+2*new_b(j)*eta(i,j))*mu;   
            end

            if rank(a_de) < size(a_de,2)
               
                break 
            end

            new_a(iind,j)=inv(a_de)*a_nu/2;  
            clear a_nu a_de
        end               
        clear j i

        if norm(new_a(1:end)-par_a(1:end),2)+norm(new_b-par_b,2)...
                +norm(Sigma(1:end)-par_Sigma(1:end),2) < 0.0001 
            converge = 1;
        end    
        n = n + 1;

    end   


    if rank(Sigma) < domain
        rsigma = par_Sigma; 
    else
        rsigma = Sigma; 
    end

    ra=new_a; rb=new_b;reta = eta; reps=xi; rsigma = Sigma; 
    mu_i = MU; sig_i = SIGMA;