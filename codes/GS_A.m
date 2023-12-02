function [U_half,V_half,P]=GS_A(U,V,P,F_U,F_V,F_D,N)%实现DGS迭代
U_A=zeros(N+1,N);
P_U=zeros(N+1,N);
U_half=zeros(N+1,N);
V_A=zeros(N,N+1);
P_V=zeros(N,N+1);
V_half=zeros(N,N+1);
h=1/N;

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
    U_temp=F_U-U_A-P_U;
    V_temp=F_V-V_A-P_V;
    
    i=1;j=1;
    U_half(i,j)=U(i,j)+U_temp(i,j);
    for i=2:N
        U_half(i,1)=1/3*(U_half(i-1,1)-U(i-1,1)+3*U(i,1)+h^2*U_temp(i,1));
    end
    i=N+1;j=1;
    U_half(i,j)=U(i,j)+U_temp(i,j);
    
    for j=2:N-1
        i=1;
        U_half(i,j)=U(i,j)+U_temp(i,j);
        for i=2:N
            U_half(i,j)=1/4*(U_half(i,j-1)+U_half(i-1,j)-(U(i,j-1)+U(i-1,j)-4*U(i,j))+h^2*U_temp(i,j));
        end
        i=N+1;
        U_half(i,j)=U(i,j)+U_temp(i,j);
    end
    
    i=1;j=N;
    U_half(i,j)=U(i,j)+U_temp(i,j);
    for i=2:N
        U_half(i,N)=1/3*(U_half(i-1,N)+U_half(i,N-1)-(U(i-1,N)-3*U(i,N)+U(i,N-1))+h^2*U_temp(i,N));
    end
    i=N+1;j=N;
    U_half(i,j)=U(i,j)+U_temp(i,j);
    
    i=1;j=1;
    V_half(i,j)=V(i,j)+V_temp(i,j);
    for j=2:N
        V_half(1,j)=1/3*(V_half(1,j-1)-(-3*V(1,j)+V(1,j-1))+h^2*V_temp(1,j));
    end
    i=1;j=N+1;
    V_half(i,j)=V(i,j)+V_temp(i,j);
    
    for i=2:N-1
        j=1;
        V_half(i,j)=V(i,j)+V_temp(i,j);
        for j=2:N
            V_half(i,j)=1/4*(V_half(i,j-1)+V_half(i-1,j)-(V(i,j-1)+V(i-1,j)-4*V(i,j))+h^2*V_temp(i,j));
        end
        j=N+1;
        V_half(i,j)=V(i,j)+V_temp(i,j);
    end
    
    i=N;j=1;
    V_half(i,j)=V(i,j)+V_temp(i,j);
    for j=2:N
        V_half(N,j)=1/3*(V_half(N-1,j)+V_half(N,j-1)-(V(N-1,j)+V(N,j-1)-3*V(N,j))+h^2*V_temp(N,j));
    end
    i=N;j=N+1;
    V_half(i,j)=V(i,j)+V_temp(i,j);
    
    for i=2:N-1
        for j=2:N-1
            r_temp=-1/h*(U_half(i+1,j)-U_half(i,j)+V_half(i,j+1)-V_half(i,j))-F_D(i,j);
            delta=r_temp*h/4;
            U_half(i,j)=U_half(i,j)-delta;
            U_half(i+1,j)=U_half(i+1,j)+delta;
            V_half(i,j)=V_half(i,j)-delta;
            V_half(i,j+1)=V_half(i,j+1)+delta;
            P(i,j)=P(i,j)+r_temp;
            P(i+1,j)=P(i+1,j)-1/4*r_temp;
            P(i-1,j)=P(i-1,j)-1/4*r_temp;
            P(i,j+1)=P(i,j+1)-1/4*r_temp;
            P(i,j-1)=P(i,j-1)-1/4*r_temp;
        end
    end
    for i=2:N-1
        j=N;
        r_temp=-1/h*(U_half(i+1,j)-U_half(i,j)+V_half(i,j+1)-V_half(i,j))-F_D(i,j);
        delta=r_temp*h/3;
        U_half(i,j)=U_half(i,j)-delta;
        U_half(i+1,j)=U_half(i+1,j)+delta;
        V_half(i,j)=V_half(i,j)-delta;
        P(i,j)=P(i,j)+r_temp;
        P(i+1,j)=P(i+1,j)-1/3*r_temp;
        P(i-1,j)=P(i-1,j)-1/3*r_temp;
        P(i,j-1)=P(i,j-1)-1/3*r_temp;
        
        j=1;
        r_temp=-1/h*(U_half(i+1,j)-U_half(i,j)+V_half(i,j+1)-V_half(i,j))-F_D(i,j);
        delta=r_temp*h/3;
        U_half(i,j)=U_half(i,j)-delta;
        U_half(i+1,j)=U_half(i+1,j)+delta;
        V_half(i,j+1)=V_half(i,j+1)+delta;
        P(i,j)=P(i,j)+r_temp;
        P(i+1,j)=P(i+1,j)-1/3*r_temp;
        P(i-1,j)=P(i-1,j)-1/3*r_temp;
        P(i,j+1)=P(i,j+1)-1/3*r_temp;
    end
    for j=2:N-1
        i=N;
        r_temp=-1/h*(U_half(i+1,j)-U_half(i,j)+V_half(i,j+1)-V_half(i,j))-F_D(i,j);
        delta=r_temp*h/3;
        V_half(i,j)=V_half(i,j)-delta;
        V_half(i,j+1)=V_half(i,j+1)+delta;
        U_half(i,j)=U_half(i,j)-delta;
        P(i,j)=P(i,j)+r_temp;
        P(i,j+1)=P(i,j+1)-1/3*r_temp;
        P(i-1,j)=P(i-1,j)-1/3*r_temp;
        P(i,j-1)=P(i,j-1)-1/3*r_temp;
        
        i=1;
        r_temp=-1/h*(U_half(i+1,j)-U_half(i,j)+V_half(i,j+1)-V_half(i,j))-F_D(i,j);
        delta=r_temp*h/3;
        V_half(i,j)=V_half(i,j)-delta;
        V_half(i,j+1)=V_half(i,j+1)+delta;
        U_half(i+1,j)=U_half(i+1,j)+delta;
        P(i,j)=P(i,j)+r_temp;
        P(i,j+1)=P(i,j+1)-1/3*r_temp;
        P(i+1,j)=P(i+1,j)-1/3*r_temp;
        P(i,j-1)=P(i,j-1)-1/3*r_temp;
    end
    
    i=1;j=1;
    r_temp=-1/h*(U_half(i+1,j)-U_half(i,j)+V_half(i,j+1)-V_half(i,j))-F_D(i,j);
    delta=r_temp*h/2;
    
    V_half(i,j+1)=V_half(i,j+1)+delta;
    U_half(i+1,j)=U_half(i+1,j)+delta;
    P(i,j)=P(i,j)+r_temp;
    P(i,j+1)=P(i,j+1)-1/2*r_temp;
    P(i+1,j)=P(i+1,j)-1/2*r_temp;
    
    i=N;j=1;
    r_temp=-1/h*(U_half(i+1,j)-U_half(i,j)+V_half(i,j+1)-V_half(i,j))-F_D(i,j);
    delta=r_temp*h/2;
    
    V_half(i,j+1)=V_half(i,j+1)+delta;
    U_half(i,j)=U_half(i,j)-delta;
    P(i,j)=P(i,j)+r_temp;
    P(i-1,j)=P(i-1,j)-1/2*r_temp;
    P(i,j+1)=P(i,j+1)-1/2*r_temp;
    
    i=1;j=N;
    r_temp=-1/h*(U_half(i+1,j)-U_half(i,j)+V_half(i,j+1)-V_half(i,j))-F_D(i,j);
    delta=r_temp*h/2;
    
    V_half(i,j)=V_half(i,j)-delta;
    U_half(i+1,j)=U_half(i+1,j)+delta;
    P(i,j)=P(i,j)+r_temp;
    P(i+1,j)=P(i+1,j)-1/2*r_temp;
    P(i,j-1)=P(i,j-1)-1/2*r_temp;
    
    i=N;j=N;
    r_temp=-1/h*(U_half(i+1,j)-U_half(i,j)+V_half(i,j+1)-V_half(i,j))-F_D(i,j);
    delta=r_temp*h/2;
    
    V_half(i,j)=V_half(i,j)-delta;
    U_half(i,j)=U_half(i,j)-delta;
    P(i,j)=P(i,j)+r_temp;
    P(i-1,j)=P(i-1,j)-1/2*r_temp;
    P(i,j-1)=P(i,j-1)-1/2*r_temp;
end




