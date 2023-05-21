%initialization
function [a0,b0,eta0,eps0,Sigma]=init(u,domain,indic)
[person,item]=size(u);

%%% Initial values.
r=identify(u);
r(find(r>0.9))=0.9;
r(find(r<0))=abs(r(find(r<0,1)));
a0=(ones(domain,1)*(r./sqrt(1-r.^2))).*indic;
b0=-norminv(sum(u)./person,0,1)./r;
Sigma = eye(domain);
theta=normrnd(0,1,person,domain);
xi=ones(person,1)*b0-theta*a0;
eta0=(exp(xi)./(1+exp(xi))-0.5)./(2*xi);    
eps0=xi;
 






