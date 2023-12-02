function [U,V,P,t]=I_UIM(F_U,F_V,N,alpha,v1,v2,tau,bottom,K_max,epsilon,error)%Uzawa迭代
r0=sqrt(norm(F_U,2)^2+norm(F_V,2)^2);
P=zeros(N,N);
U=zeros(N+1,N);
V=zeros(N,N+1);
t=0;
while 1
    t=t+1;
    P_temp=apply_BT(U,V,N);
    error2=tau*norm(P_temp);
    [Utemp,Vtemp]=apply_B(P,N);
    F_U1=F_U-Utemp;F_V1=F_V-Vtemp;
    [U,V,s,v]=v_cycle_PCG(F_U1,F_V1,N,K_max,epsilon,error,bottom,v1,v2,error2);
    fprintf("PCG ite=%d \n",s);
    disp(v);
    P1=apply_BT(U,V,N);
    P=P+alpha*P1;
    [U_A,V_A,P_A]=apply_matrix(U,V,P,N);
    if sqrt(norm(F_U-U_A,2)^2+norm(F_V-V_A,2)^2++norm(P_A,2)^2)/r0<1e-8
        break
    end
end