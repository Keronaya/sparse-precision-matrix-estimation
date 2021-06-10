function Sigma=yangbenxiefancha(A)
for i=1:size(A,2) 
    for j=1:size(A,2) 
        Sigma(i,j)=sum(A(:,i).*A(:,j))/size(A,1);
    end
end