function [X iter]=shrink_G(A,rho,positive)
D=size(A,1);
X=sign(A).*max(abs(A)-rho+rho*eye(D),0);
iter=0;
if positive==1 & min(eig(X))<0

X0=X+10^-5;
Lambda=zeros(D,D);
mu=1;

while norm(X0-X,'fro')>10^-6
[Utemp,Stemp]=eig(X0+mu*Lambda);
X=Utemp*max(Stemp,0)*Utemp';
temp=mu*(A-Lambda)+X;
X0=(sign(temp).*max(abs(temp)-rho*mu+rho*mu*eye(D),0))/(1+mu);
Lambda=Lambda-(X-X0)/mu;
iter=iter+1;
end
end