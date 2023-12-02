function [U,V,s]=V_cycle(F_U,F_V,v1,v2,L,bottom)%通过V-cycle多重网格方法，配合DGS磨光子求解
K=log(L/bottom)/log(2);
N=L;
F_P=zeros(N,N);
U_store=cell(1,K+1);
V_store=cell(1,K+1);
P_store=cell(1,K+1);
F_U_store=cell(1,K);
F_V_store=cell(1,K);
F_P_store=cell(1,K);

U_store{1}=zeros(N+1,N);
V_store{1}=zeros(N,N+1);
P_store{1}=zeros(N,N);
F_U_store{1}=F_U;
F_V_store{1}=F_V;
F_P_store{1}=F_P;
r0=sqrt(norm(F_U,2)^2+norm(F_V,2)^2);
s=0;
while 1
    s=s+1;
    for i=1:K+1
        if i~=1
            U=zeros(N+1,N);
            V=zeros(N,N+1);
            P=zeros(N,N);
        else
            U=U_store{1};
            V=V_store{1};
            P=P_store{1};
            F_U=F_U_store{1};
            F_V=F_V_store{1};
            F_P=F_P_store{1};
        end
        if i~=K+1
            for k=1:v1
                [U,V,P]=GS_A(U,V,P,F_U,F_V,F_P,N);
            end
        else
            for k=1:50
                [U,V,P]=GS_A(U,V,P,F_U,F_V,F_P,N);
            end
            
        end
        U_store{i}=U;
        V_store{i}=V;
        P_store{i}=P;
        if i~=1 && i~=K+1
            F_U_store{i}=F_U;
            F_V_store{i}=F_V;
            F_P_store{i}=F_P;
        end
        if i~=K+1
            [U,V,P]=apply_matrix(U,V,P,N);
            R_U=F_U-U;R_V=F_V-V;R_P=F_P-P;
            [F_U,F_V,F_P]=limit(R_U,R_V,R_P,N);
            N=N/2;
        end
    end
    %N为底层的大小
    for j=1:K
        i=K+1-j;
        [U,V,P]=lift(U_store{i+1},V_store{i+1},P_store{i+1},N);
        N=N*2;
        U=U_store{i}+U;
        V=V_store{i}+V;
        P=P_store{i}+P;
        for k=1:v2
            [U,V,P]=GS_A(U,V,P,F_U_store{i},F_V_store{i},F_P_store{i},N);
        end
        U_store{i}=U;
        V_store{i}=V;
        P_store{i}=P;
    end
    [U1,V1,P1]=apply_matrix(U,V,P,N);
    rh=sqrt(norm(U1-F_U_store{1})^2+norm(V1-F_V_store{1})^2+norm(P1,2)^2);
    if rh/r0<1e-8
        break
    end
end
end