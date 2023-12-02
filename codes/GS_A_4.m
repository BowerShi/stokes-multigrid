function [U,V]=GS_A_4(U,V,F_U,F_V,N)%关于矩阵A的对称GS迭代
h=1/N;
for i=2:N
    U(i,1)=1/3*(U(i-1,1)+U(i+1,1)+U(i,2)+h^2*F_U(i,1));
end
for j=2:N-1
    for i=2:N
        U(i,j)=1/4*(U(i,j+1)+U(i,j-1)+U(i-1,j)+U(i+1,j)+h^2*F_U(i,j));
    end
end
for i=2:N
    U(i,N)=1/3*(U(i-1,N)+U(i+1,N)+U(i,N-1)+h^2*F_U(i,N));
end

for j=2:N
    V(1,j)=1/3*(V(2,j)+V(1,j+1)+V(1,j-1)+h^2*F_V(1,j));
end
for i=2:N-1
    for j=2:N
        V(i,j)=1/4*(V(i,j+1)+V(i,j-1)+V(i-1,j)+V(i+1,j)+h^2*F_V(i,j));
    end
end
for j=2:N
    V(N,j)=1/3*(V(N-1,j)+V(N,j+1)+V(N,j-1)+h^2*F_V(N,j));
end

for i=N:-1:2
    U(i,N)=1/3*(U(i-1,N)+U(i+1,N)+U(i,N-1)+h^2*F_U(i,N));
end

for j=N-1:-1:2
    for i=N:-1:2
        U(i,j)=1/4*(U(i,j+1)+U(i,j-1)+U(i-1,j)+U(i+1,j)+h^2*F_U(i,j));
    end
end
for i=N:-1:2
    U(i,1)=1/3*(U(i-1,1)+U(i+1,1)+U(i,2)+h^2*F_U(i,1));
end

for j=N:-1:2
    V(N,j)=1/3*(V(N-1,j)+V(N,j+1)+V(N,j-1)+h^2*F_V(N,j));
end

for i=N-1:-1:2
    for j=N:-1:2
        V(i,j)=1/4*(V(i,j+1)+V(i,j-1)+V(i-1,j)+V(i+1,j)+h^2*F_V(i,j));
    end
end

for j=N:-1:2
    V(1,j)=1/3*(V(2,j)+V(1,j+1)+V(1,j-1)+h^2*F_V(1,j));
end



