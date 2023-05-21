function r=identify(u)
    [person,item]=size(u);
    scale=sum(u,2);  
    p = mean(u,1);   % proportion of positive responses for each item
    q=1-p; % proportion of negative responses for each item
    y=norminv(p,0,1); % returns the inverse of the standard normal cumulative distribution function (cdf), evaluated at the probability values in p.
    y=normpdf(y,0,1); % 
    s=std(scale);
    for i=1:item
        x1=scale(u(:,i)==1);  % num of examinee who has a correct response for item i
        x2=scale(u(:,i)==0); % num of examinee who has a incorrect response for item i
        r(i)=(mean(x1)-mean(x2))/s*p(i)*q(i)/y(i);
    end 