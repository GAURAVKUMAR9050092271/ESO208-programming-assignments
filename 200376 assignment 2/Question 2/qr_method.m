function qr_method(A,itr,tol)
[n,~] = size(A);
Q = zeros(n,n);
e = zeros(n,1);
R = zeros(n,n);
A_prev=zeros(n,1);
err= 100;
t=0;
while (err>tol && t<itr)
    for p=1:1:n
        Q(p,1) = A(p,1)/norm(A(:,1));
    end
    for j=2:n
        sum=0;
        for i=1:j-1
            sum =sum + (Q(:,i)'*A(:,j))*Q(:,i);
        end
        e(:,1) = A(:,j)- sum;
        for i=1:n
            Q(i,j) = e(i,1)/norm(e);
        end
    end
    for i=1:n
        for j=1:n
            R(i,j) = Q(:,i)'*A(:,j);
        end
    end
    A = R*Q;
    er = zeros(n,1);
    for i=1:n
        er(i)= abs((A(i,i)-A_prev(i))/A(i,i))*100;
        A_prev(i)=A(i,i);
    end
    err=max(er);
    t=t+1;
end
B = zeros(n,1);
for i=1:1:n
    B(i,1)= A(i,i);
end
prt = fopen("qr_out.txt",'w');
fprintf(prt,"The eigen values of given matrix is: \n");
fprintf(prt,"%0.4f\n",B);
fprintf(prt,"Iterations: %d\n",t);
fclose(prt);

end


