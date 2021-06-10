a=10;
for i =1:1:a
A=0;
M=0;
N=50;p=100;n=100;para=0.3;
B=zeros(p,p);
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
 B=Theta+B;
end
B=B/a;
for i=1:1:p
    for j=1:1;p
        if(Theta0(i,j)~=0)
            A=(B(i,j)-Theta0(i,j))^2+A;
        else
            M=(B(i,j)-Theta0(i,j))^2+M;
        end
    end
end
D=A/(3*p-2)
M=M/(p^2-(3*p-2))
            
