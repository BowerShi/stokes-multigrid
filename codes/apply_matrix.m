function [U1,V1,P1]=apply_matrix(U,V,P,N)
h=1/N;
U_A=zeros(N+1,N);
P_U=zeros(N+1,N);

V_A=zeros(N,N+1);
P_V=zeros(N,N+1);

P1=zeros(N,N);

for j=2:N-1
    for i=2:N
        U_A(i,j)=-1/(h^2)*(U(i,j+1)+U(i,j-1)+U(i-1,j)+U(i+1,j)-4*U(i,j));
        P_U(i,j)=1/h*(P(i,j)-P(i-1,j));
    end
end
for i=2:N
    U_A(i,1)=-1/(h^2)*(U(i-1,1)+U(i+1,1)-2*U(i,1)+U(i,2)-U(i,1));
    P_U(i,1)=1/h*(P(i,1)-P(i-1,1));
    U_A(i,N)=-1/(h^2)*(U(i-1,N)+U(i+1,N)-2*U(i,N)-U(i,N)+U(i,N-1));
    P_U(i,N)=1/h*(P(i,N)-P(i-1,N));
end
for j=1:N
    U_A(1,j)=U(1,j);
    U_A(N+1,j)=U(N+1,j);
end

for i=2:N-1
    for j=2:N
        V_A(i,j)=-1/(h^2)*(V(i,j+1)+V(i,j-1)+V(i-1,j)+V(i+1,j)-4*V(i,j));
        P_V(i,j)=1/h*(P(i,j)-P(i,j-1));
    end
end
for j=2:N
    V_A(1,j)=-1/(h^2)*(V(2,j)-V(1,j)+V(1,j+1)+V(1,j-1)-2*V(1,j));
    P_V(1,j)=1/h*(P(1,j)-P(1,j-1));
    V_A(N,j)=-1/(h^2)*(V(N-1,j)-V(N,j)+V(N,j+1)+V(N,j-1)-2*V(N,j));
    P_V(N,j)=1/h*(P(N,j)-P(N,j-1));
end
for i=1:N
    V_A(i,1)=V(i,1);
    V_A(i,N+1)=V(i,N+1);
end
for i=1:N
    for j=1:N
        P1(i,j)=-1/h*(U(i+1,j)-U(i,j)+V(i,j+1)-V(i,j));
    end
end
U1=U_A+P_U;
V1=V_A+P_V;
end