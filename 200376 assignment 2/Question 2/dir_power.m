function dir_power(A,itr,tol)
n = length(A);
x = ones(n,1);
eig = 0;
j = 0;
err= 100;
 while (j<=itr && err>=tol) 
   y=A*x; 
   s=max(abs(y));
   x=y/s;
   err=abs(eig-s)/s*100;
   eig=s;
   j = j+1;
 end
p = 1/norm(x);
x = x*p;
prt = fopen("output_direct.txt",'w');
fprintf(prt,"The maximum eigen value of given matrix is: %0.4f\n",eig);
fprintf(prt,"The corresponding eigen vector is: \n");
fprintf(prt,"%0.4f\n",x);
fprintf(prt,"Iterations: %d\n",j);
fclose(prt);
end
    