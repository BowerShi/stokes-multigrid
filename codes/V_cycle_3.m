function [U,V,s]=V_cycle_3(F_U,F_V,v1,v2,error,L,bottom)%多重网格方法配合对称GS迭代求解问题3
K=log(L/bottom)/log(2);
N=L;
U_store=cell(1,K+1);
V_store=cell(1,K+1);
F_U_store=cell(1,K);
F_V_store=cell(1,K);
U_store{1}=zeros(N+1,N);
V_store{1}=zeros(N,N+1);
F_U_store{1}=F_U;
F_V_store{1}=F_V;
r0=sqrt(norm(F_U_store{1})^2+norm(F_V_store{1})^2);
s=0;
ts=0;
while 1
    s=s+1;
    for i=1:K+1
        if i~=1
            U=zeros(N+1,N);
            V=zeros(N,N+1);
        else
            U=U_store{1};
            V=V_store{1};

            F_U=F_U_store{1};
            F_V=F_V_store{1};
        end
        if i~=K+1
            t1=toc;
            for t=1:v1
                [U,V]=GS_A_4(U,V,F_U,F_V,N);
            end
            t2=toc;
            ts=ts+t2-t1;
        else
            [U,V]=cal_A_3(F_U,F_V,N);
        end
        U_store{i}=U;
        V_store{i}=V;
        if i~=1 && i~=K+1
            F_U_store{i}=F_U;
            F_V_store{i}=F_V;
        end
        if i~=K+1
            [U,V]=apply_A(U,V,N);
            R_U=F_U-U;R_V=F_V-V;
            [F_U,F_V]=limit_3(R_U,R_V,N);
            N=N/2;
        end
    end
    %N为底层的大小
    for j=1:K
        i=K+1-j;
        [U,V]=lift_3(U_store{i+1},V_store{i+1},N);
        N=N*2;
        U=U_store{i}+U;
        V=V_store{i}+V;
        t1=toc;
        for t=1:v2
            [U,V]=GS_A_4(U,V,F_U_store{i},F_V_store{i},N);
        end
        t2=toc;
        ts=ts+t2-t1;
        U_store{i}=U;
        V_store{i}=V;
    end
    [U1,V1]=apply_A(U,V,N);
    rh=sqrt(norm(U1-F_U_store{1})^2+norm(V1-F_V_store{1})^2);
    if rh/r0<error
        break
    end
end
end