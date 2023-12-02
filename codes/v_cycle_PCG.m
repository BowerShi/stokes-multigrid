function [U,V,s,v]=v_cycle_PCG(F_U,F_V,N,K_max,epsilon,error,bottom,v1,v2,error2)%V-cycle作为预条件子的预优共轭梯度法
U=0*F_U;
V=0*F_V;
s=0;R_U=F_U;R_V=F_V;
temp1=sqrt(norm(F_U,2)^2+norm(F_V,2)^2);
v=zeros(2,1);
while sqrt(norm(R_U,2)^2+norm(R_V,2)^2)>max(epsilon*temp1,error2) && s<K_max
    time1=toc;
    [U_z,V_z,v_cycle_num]=V_cycle_3(R_U,R_V,v1,v2,error,N,bottom);
    time2=toc;
    s=s+1;
    v(1,s)=v_cycle_num;
    v(2,s)=time2-time1;
    if s==1
        P_U=U_z; P_V=V_z;
        temp2=v_cycle_matrix_product(R_U,U_z);
        temp3=v_cycle_matrix_product(R_V,V_z);
        rho=temp2+temp3;
    else
        rho_tilda=rho;
        temp2=v_cycle_matrix_product(R_U,U_z);
        temp3=v_cycle_matrix_product(R_V,V_z);
        rho=temp2+temp3;
        beta=rho/rho_tilda; P_U=U_z+P_U*beta; P_V=V_z+P_V*beta;
    end
    [U_w,V_w]=apply_A(P_U,P_V,N);
    temp2=v_cycle_matrix_product(P_U,U_w);
    temp3=v_cycle_matrix_product(P_V,V_w);
    temp4=temp2+temp3;
    alpha=rho/temp4;
    U=U+alpha*P_U; V=V+alpha*P_V;
    R_U=R_U-alpha*U_w; R_V=R_V-alpha*V_w;
end
end