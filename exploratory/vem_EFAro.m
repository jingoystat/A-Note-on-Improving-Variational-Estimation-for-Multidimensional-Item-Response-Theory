%GVEM algorithm for 2PL model - EFA. 
function [ra, rb, rsigma, reps, reta, mu_i, sig_i, n, converge]=vem_EFAro(u,domain,new_a,new_b,eta)
[person,item]=size(u);
converge=1;
xi = eta;  
n = 0;

while converge == 1 
        for i = 1:person 
        sigma_part=zeros(domain,domain);mu_part=zeros(domain,1);
        for j=1:item
            sigma_part=sigma_part+eta(i,j)*new_a(:,j)*new_a(:,j)';   
            mu_part=mu_part+(2*eta(i,j)*new_b(j)+u(i,j)-0.5)*new_a(:,j);     
        end
        sigmahat=inv(eye(domain)+2*sigma_part);    
        muhat=sigmahat*mu_part; 
        SIGMA(:,:,i)=sigmahat;
        MU(:,i)=muhat;
        xi(i,:)=sqrt(new_b.^2-2*new_b.*(new_a'*muhat)'+diag(new_a'*(sigmahat+muhat*muhat')*new_a)');
    end
    
    eta=(exp(xi)./(1+exp(xi))-0.5)./(2*xi);  
    for i=1:person
        for j=1:item
            if abs(xi(i,j)) < 0.01
                eta(i,j) = 0.125;  
            end
        end    
    end
    clear i j 
    
    
    par_b=new_b;
    b_part = zeros(person,item);
    for i=1:person
        b_part(i,:)=(new_a'*MU(:,i))';  
    end
    new_b=sum(0.5-u+2*eta.*b_part)./sum(2*eta);  
    clear i
    
   
    par_a=new_a;
    for j=1:item
        a_nu=zeros(domain,1); a_de=zeros(domain,domain);
        for i=1:person
            sigma=SIGMA(:,:,i);
            a_de=a_de+eta(i,j)*sigma+eta(i,j)*MU(:,i)*MU(:,i)';
            a_nu=a_nu+(u(i,j)-0.5+2*new_b(j)*eta(i,j))*MU(:,i);   
        end
        new_a(:,j)=inv(a_de) * a_nu/2;  
        clear a_nu a_de
    end
    clear j i
        
    if norm(new_b-par_b,2) + norm(new_a(1:end)-par_a(1:end),2) <0.0001  
        converge=0;
    end
    n = n+1;
end



[B,T]= rotatefactors(new_a','method','promax');
new_a = B'; 
Sigma = inv(T'*T);
mu_i = inv(T) * MU;
for i = 1:person
   sig_i(:,:,i) = inv(T) * SIGMA(:,:,i) * inv(T)';
end
ra=new_a; rb=new_b; reta = eta; reps=xi; rsigma = Sigma;

