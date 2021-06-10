for i =1:1:5
B=0;
N=50;p=100;n=100;para=0.3;
C=zeros(p*p,N);
for i=1:N
    Theta0= tridiag(p,1,para,para);
    lambda=sqrt(log(p)/n);
    rho=1;
    positive=1;
    u=zeros(p,1)';
    R=inv(Theta0);
    Y=multivrandn(u,R,n);
    Sigma=yangbenxiefancha(Y);
    [X Y iter positivity_fail]=inverse_cov(Sigma,lambda,rho,positive);
    Theta=X;
  
end
for i=1:1:p
    for j=1:1:p
        A=(sum((Theta(i,j)-Theta0(i,j))^2))/(p^2);
    end 
end
B=B+A;
end
D=B/N