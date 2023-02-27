fprintf("Enter the method you want to use\n");
method=input('1. Least square fitting polynomial\n2. Lagrange polynomial\n3. Cubic spline\n','s');
if method=='1'
    least()
else
    if method=='2'
        lagrange()
    else
        cubicspline()
    end
end
    