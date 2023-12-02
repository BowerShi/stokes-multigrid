function P1=apply_BT(U,V,N)
h=1/N;
P1=zeros(N,N);
for i=1:N
    for j=1:N
        P1(i,j)=-1/h*(U(i+1,j)-U(i,j)+V(i,j+1)-V(i,j));
    end
end
end