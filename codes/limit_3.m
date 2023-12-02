function [U1,V1]=limit_3(U,V,N)
U1=zeros(N/2+1,N/2);
V1=zeros(N/2,N/2+1);
for j=1:N/2
    U1(1,j)=1/2*(U(1,2*j)+U(1,2*j-1));
    U1(N/2+1,j)=1/2*(U(N+1,2*j)+U(N+1,2*j-1));
end
for i=2:N/2
    for j=1:N/2
        U1(i,j)=1/4*(U(2*i-1,2*j)+U(2*i-1,2*j-1))+1/8*(U(2*i-2,2*j)+U(2*i,2*j)+U(2*i-2,2*j-1)+U(2*i,2*j-1));
    end
end
for j=1:N/2
    V1(j,1)=1/2*(V(2*j,1)+V(2*j-1,1));
    V1(j,N/2+1)=1/2*(V(2*j,N+1)+V(2*j-1,N+1));
end
for j=2:N/2
    for i=1:N/2
        V1(i,j)=1/4*(V(2*i,2*j-1)+V(2*i-1,2*j-1))+1/8*(V(2*i,2*j-2)+V(2*i-1,2*j-2)+V(2*i,2*j)+V(2*i-1,2*j));
    end
end
end