function [a,b,u,theta,Sigma]=simulation_2PL(person, item, domain, indic, sigL, sigR, model)
    % This function is used to generate response data for MIRT models.
    % indic is the indicator of multidimensional structure.
    % The size of a is domain*item, the size of b is 1*item,the size of theta is person*(domain+1).
    % sigL, sigR : 0.1 (low correlation), 0.7 (high correlation)


    d=1;
    P=@(theta,a,b) 1./(1+exp(-d*(theta*a)+d*kron(ones(size(theta,1),1),b)));  %person*item

    % ITEM PARAMETERS GENERATION
    %if strcmp(model,'within')
    if model == 1
        % within item model
        a=unifrnd(1,2,domain,item).*indic; %lognrnd(0,0.5,domain,item).*indic;   
        b=normrnd(0,1,1,item);
    else
        % between item model 
        a=linspace(1,2,domain)'*ones(1,item).*indic;   % domain*item %CFA
        b=normrnd(0,1,1,item);
    end

    q=1; % Index to judge whether Sigma is a positive matrix.
    while q~=0
        sig_temp=unifrnd(sigL,sigR,domain*(domain-1)/2,1);
        temp(1)=1; temp(2:domain)=[domain-1:-1:1];
        for k=1:domain-1
            Sigma(k,k)=1;
            Sigma(k+1:domain,k)=sig_temp(sum(temp(1:k)):sum(temp(1:k+1))-1,1);
            Sigma(k,k+1:domain)=sig_temp(sum(temp(1:k)):sum(temp(1:k+1))-1,1)';
        end
        Sigma(domain,domain)=1;
        [r,q]=chol(Sigma);
    end
    clear r

    theta=mvnrnd(zeros(person,domain),Sigma); % person * domain
    p=P(theta,a,b);
    u=binornd(1,p);
