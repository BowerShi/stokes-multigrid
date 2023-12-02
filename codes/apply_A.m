function [U_A,V_A]=apply_A(U,V,N)
h=1/N;
U_A=zeros(N+1,N);
V_A=zeros(N,N+1);

for j=2:N-1
    for i=2:N
        U_A(i,j)=-1/(h^2)*(U(i,j+1)+U(i,j-1)+U(i-1,j)+U(i+1,j)-4*U(i,j));
    end
end
for i=2:N
    U_A(i,1)=-1/(h^2)*(U(i-1,1)+U(i+1,1)-2*U(i,1)+U(i,2)-U(i,1));
    U_A(i,N)=-1/(h^2)*(U(i-1,N)+U(i+1,N)-2*U(i,N)-U(i,N)+U(i,N-1));
end
for j=1:N
    U_A(1,j)=U(1,j);
    U_A(N+1,j)=U(N+1,j);
end

for i=2:N-1
    for j=2:N
        V_A(i,j)=-1/(h^2)*(V(i,j+1)+V(i,j-1)+V(i-1,j)+V(i+1,j)-4*V(i,j));
    end
end
for j=2:N
    V_A(1,j)=-1/(h^2)*(V(2,j)-V(1,j)+V(1,j+1)+V(1,j-1)-2*V(1,j));
    V_A(N,j)=-1/(h^2)*(V(N-1,j)-V(N,j)+V(N,j+1)+V(N,j-1)-2*V(N,j));
end
for i=1:N
    V_A(i,1)=V(i,1);
    V_A(i,N+1)=V(i,N+1);
end
end