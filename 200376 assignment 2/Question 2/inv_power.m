function inv_power(p,itr,tol)
n = length(p);
x = ones(n,1);
eig = 0;
j = 0;
err= 100;
A = inv(p);
while (j<itr && err>=tol) 
   y=A*x; 
   s=max(abs(y));
   x=y/s;
   err=abs(eig-s)/s*100;
   eig=s;
   j = j+1;
end
 eig = 1/eig;
p = 1/norm(x);
x = x*p;
prt = fopen("inverse_out.txt",'w');
fprintf(prt,"The minimmum eigen value of given matrix is: %0.4f\n",eig);
fprintf(prt,"The corresponding eigen vector is: \n");
fprintf(prt,"%0.4f\n",x);
fprintf(prt,"Iterations: %d\n",j);
fclose(prt);
end
    
