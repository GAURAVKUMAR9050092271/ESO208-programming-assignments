function user_interface()
clc
fID = fopen("input.txt",'rt');
sizen = 1;
n = fscanf(fID,"%d",sizen);
sizeA = [n,n];
 p= fscanf(fID,"%f",sizeA);
sizei=1;
itr = fscanf(fID,"%f",sizei);
tol = fscanf(fID,"%f",sizen);
shift = fscanf(fID,"%f",sizen);
fclose(fID);
fprintf("1:Direct Power Method\n2:Inverse Power Method\n3:Shifted Power Method\n4:QR Method\n");
ch = input("Enter Choice: ");
if ch==1
    dir_power(p,itr,tol);
elseif ch==2
    inv_power(p,itr,tol);
elseif ch==3
    shift_power(p,itr,tol,shift);
elseif ch==4
    qr_method(p,itr,tol);
end
end