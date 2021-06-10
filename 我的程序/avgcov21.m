N=100;p=100;n=100;para=0.2;
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
    C(:,i)=zhixinqvjian2(Sigma,Theta,Theta0,n,0.05);
end
c1=0;c2=0;
for i=1:p*p
    for j=1:N
        if Theta0(i)~=0&C(i,j)==1
            c1=c1+1;
        elseif Theta0(i)==0&C(i,j)==1
            c2=c2+1;
        else
            c1=c1+0;
            c2=c2+0;
        end
    end
end
s1=c1/((3*p-2)*N)
s2=c2/((p*p-3*p+2)*N)
