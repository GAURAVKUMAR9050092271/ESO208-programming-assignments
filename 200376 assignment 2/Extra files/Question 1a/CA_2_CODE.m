
n=input("INPUT NO OF UNKNOWN VARIABLES");
x=input("INPUT 1 FOR GE(without pivoting)\n INPUT 2 FOR GE(with pivoting)\n INPUT 3 FOR GE(with s&p)\n INPUT 4 FOR LU BY GE(without pivoting) \n INPUT 5 FOR LU BY GE(with pivoting) \n INPUT 6 FOR LU BY CROUT \n INPUT 7 FOR CHOLESKY");
A=input("INPUT AUGUMENTED MATRIX");
if(x==1)
     l = zeros(n); 
    for k=1:n
        for i=k+1:n
                l(i,k)=A(i,k)/A(k,k);
                A(i,:)=A(i,:)- l(i,k)*A(k,:);   
        end
    end
    disp('Final Matrix is: ')
    disp(A)
     b = A(:,n+1);
    A = A(:,(1:n));
   
  
    x=zeros(n,1);
    x(n)= b(n)/A(n,n);
    
    for i=n-1:-1:1
        sum=0;
        for j=i+1:n
            sum = sum + A(i,j)*x(j);
        end
        x(i)=(b(i)-sum)/A(i,i);
    end
    
    disp('The Unknowns are: ')
    disp(x)
end
        

if(x==2)
    l = zeros(n); 
    for k=1:n
        [M,p]=max(abs(A(k:n,k)));
        P = eye(n);
        P( [k, k+p-1], : ) = P( [k+p-1, k], : );
        disp('The permutation Matrix')
        disp(P)
        disp('The Augmented Matrix after Partial Pivoting')
        A = P*A;
        
        if (A(k,k)==0)
            disp('ERROR!!')
            return;
        end
        
        for i=k+1:n
            l(i,k)=A(i,k)/A(k,k);
            A(i,:)=A(i,:)- l(i,k)*A(k,:); 
        end
    end
    disp('The Reduced Matrix is: ')
    disp(A)
    b = A(:,n+1);
    A = A(:,(1:n));
    
    
    x=zeros(n,1);
    x(n)= b(n)/A(n,n);
    
    for i=n-1:-1:1
        sum=0;
        for j=i+1:n
            sum = sum + A(i,j)*x(j);
        end
        x(i)=(b(i)-sum)/A(i,i);
    end
    
    disp('The Unknowns are: ')
    disp(x)
    
 
    
    
end
    
if (x==3)
    for i=1:n
        s1(i,1) = max(A(i,:));
        A(i,:) = A(i,:)/s1(i,1);
    end
    x1 = zeros(n,1);   
    L = eye(n);
    for i=1:n-1
        [m,ind]=max(abs(A(i:n,i)));
        temp1=A(i,:);
        A(i,:)=A(ind+i-1,:);
        A(ind+i-1,:)=temp1;
        for j=i+1:n
            l1 = A(j,i)/A(i,i);
            A(j,:) = A(j,:) - l1.*A(i,:);
        end
    end
    disp("Final Matrix");
    disp(A);
    
     x(n) = A(n,n+1)/A(n,n);
    for i=n-1:-1:1
        sum = 0;
        for j=i+1:n
            sum = sum + A(i,j)*x(j);
        end
        x(i) = (A(i,n+1) - sum)/A(i,i);
    end
    disp('The Unknowns are: ')
    disp(x);
end
    
    
if(x==4)
   L = zeros(n); 
    for k=1:n
        L(k,k)=1;
        if (A(k,k)==0)
            disp('Pivoting  is Required!')
            return;
        end
        
        for i=k+1:n
                L(i,k)=A(i,k)/A(k,k);
                A(i,:)=A(i,:)- L(i,k)*A(k,:);   
        end
    end
    b = A(:,n+1);
    A = A(:,(1:n));    
    disp('The Upper Triangular Matrix')
    U=A;
    disp(U);
    disp('The Lower Triangular Matrix')
    disp(L)

    x=zeros(n,1);
    x(n)= b(n)/A(n,n);
    
    for i=n-1:-1:1
        sum=0;
        for j=i+1:n
            sum = sum + A(i,j)*x(j);
        end
        x(i)=(b(i)-sum)/A(i,i);
    end
    
    disp('The Unknowns are: ')
    disp(x)
end




if(x==5)
    a=A(:,1:n);
    b=A(:,n+1);
    for i=1:n-1
    [m,p]=max(abs(A(i:n,i)));
    c=A(i,:);
    A(i,:)=A(p+i-1,:);
    A(p+i-1,:)=c;
    end
    L=zeros(n);
    U=zeros(n);
    for i = 1:n
       L(i,i) = 1;
    end
    for k=1:n
        for j=k:n
            s2=0;
            for m=1:k-1
                s2=s2+L(k,m)*U(m,j);
            end
            U(k,j)=(A(k,j)-s2);
        end    
        for i=k+1:n
            s1=0;
            for m=1:k-1
                s1=s1+L(i,m)*U(m,k);
            end
            L(i,k)=(A(i,k)-s1)/U(k,k);
        end
    end  
    disp('The Upper Triangular Matrix')
    disp(U);
    disp('The Lower Triangular Matrix')
    disp(L)
    %L,U
    y=zeros(n,1);
    y(1)=b(1)/L(1,1);
    for k=2:n
     y(k)=(b(k)-L(k,1:k-1)*y(1:k-1))/L(k,k);
    end
    x=zeros(n,1);
    x(n)=y(n)/U(n,n);
    for k=n-1:-1:1
        x(k)=(y(k)-U(k,k+1:n)*x(k+1:n))/U(k,k);
    end
    disp("The unknowns x are");
    disp(x);
end



if(x==6)
    L=zeros(n);
    U=zeros(n);
    for i = 1:n
       U(i,i) = 1;
    end
    for k=1:n
        for i=k:n
            s1=0;
            for m=1:k-1
                s1=s1+L(i,m)*U(m,k);
            end
            L(i,k)=A(i,k)-s1;
        end
        for j=k+1:n
            s2=0;
            for m=1:k-1
                s2=s2+L(k,m)*U(m,j);
            end
            U(k,j)=(A(k,j)-s2)/L(k,k);
        end    
    end
    disp('The Upper Triangular Matrix');
    disp(U);
    disp('The Lower Triangular Matrix');
    disp(L);
    x = zeros(n,1);
    y = zeros(n,1);
    y(1) = A(1,n+1)/L(1,1);
    for i=2:n
        sum = 0;
        for j=1:i
            sum = sum + L(i,j)*y(j,:);
        end
        y(i,:) = (A(i,n+1) - sum)/L(i,i);
    end
    x(n) = y(n,:)/U(n,n);
    for i=n-1:-1:1
        sum = 0;
        for j=i+1:n
            sum = sum + U(i,j)*x(j,:);
        end
        x(i,:) = (y(i,:) - sum)/U(i,i);
    end
    disp(x)
end



if(x==7)
    L=zeros(n);
    for i=1:n
        
        for j=i:n
            if (j==i)
                L(j,i)=sqrt(A(i,j) - sum(L(i,1:(i-1)).*L(j,1:(i-1)) ));
            else
                L(j,i)=A(i,j) - sum(L(i,1:(i-1)).*L(j,1:(i-1)) );
            end
        end
    end
    disp('The Chloskey Factor, Lc is');
    disp(L);
    
 
end