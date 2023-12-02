function [P_U,P_V]=apply_B(P,N)
h=1/N;
P_U=zeros(N+1,N);
P_V=P_U';
for j=2:N-1
    for i=2:N
        P_U(i,j)=1/h*(P(i,j)-P(i-1,j));
    end
end
for i=2:N
    P_U(i,1)=1/h*(P(i,1)-P(i-1,1));
    P_U(i,N)=1/h*(P(i,N)-P(i-1,N));
end
for i=2:N-1
    for j=2:N
        P_V(i,j)=1/h*(P(i,j)-P(i,j-1));
    end
end
for j=2:N
    P_V(1,j)=1/h*(P(1,j)-P(1,j-1));
    P_V(N,j)=1/h*(P(N,j)-P(N,j-1));
end
end