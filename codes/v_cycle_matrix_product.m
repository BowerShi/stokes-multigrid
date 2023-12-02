function sum=v_cycle_matrix_product(A,B)%把矩阵当成向量求内积
sum=0;
[n1,n2]=size(A);
for i=1:n1
    for j=1:n2
        sum=sum+A(i,j)*B(i,j);
    end
end
end