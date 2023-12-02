function [U1,V1,P1]=lift(U,V,P,N)
U1=zeros(2*N+1,2*N);
V1=zeros(2*N,2*N+1);
P1=zeros(2*N,2*N);
for i=1:N+1
    for j=1:N
        U1(2*i-1,2*j-1)=U(i,j);
        U1(2*i-1,2*j)=U(i,j);
    end
end
for i=1:N
    for j=1:N
        U1(2*i,2*j-1)=1/2*(U(i,j)+U(i+1,j));
        U1(2*i,2*j)=1/2*(U(i,j)+U(i+1,j));
    end
end
for j=1:N+1
    for i=1:N
        V1(2*i,2*j-1)=V(i,j);
        V1(2*i-1,2*j-1)=V(i,j);
    end
end
for j=1:N
    for i=1:N
        V1(2*i,2*j)=1/2*(V(i,j)+V(i,j+1));
        V1(2*i-1,2*j)=V1(2*i,2*j);
    end
end
for j=1:N
    for i=1:N
        P1(2*i-1,2*j-1)=P(i,j);
        P1(2*i,2*j-1)=P(i,j);
        P1(2*i-1,2*j)=P(i,j);
        P1(2*i,2*j)=P(i,j);     
    end
end
end
