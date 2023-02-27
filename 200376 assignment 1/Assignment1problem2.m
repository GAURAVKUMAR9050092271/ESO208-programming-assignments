clear
close all
clc

% program starts
method = input('Choose the method of solution by selecting number shown with method:\n Muller Method-1,\n Bairstow Method-2 \n','s');

% Muller Method
if(method == '1')
    funcInput = input("Enter the polynomial whose roots are to be calculated :\n", 's');
    funcInput = strcat("@(x)", funcInput);
    func = str2func(funcInput);
    x0 = input("What is the  first starting value for initializing the algorithm:\n");
    x1 = input("What is the second starting value for initializing the algorithm:\n");
    x2 = input("What is the third starting value for initializing the algorithm:\n");
    disp("Now Enter the stopping criteria");
    m_iter = input("What is the allowed  maximum number of iterations for the algorithm:\n");
    max_error_perc = input("What is the maximum percentage error allowed:\n");

    iter = 0;
    i = 1;
    c = func(x2);
    partialdifferentials0 = (func(x1) - func(x0))/(x1 - x0);
    partialdifferentials1 = (func(x2) - func(x1))/(x2 - x1);
    partial1 = (x2 - x1);
    partial0 = (x1 - x0);
    a = (partialdifferentials1 - partialdifferentials0)/(partial0 + partial1);
    b = a*partial1 + partialdifferentials1;
    error = 100;
    discriminant = 2*c/(b + sqrt(b*b - 4*a*c));
    x3 = x2 - discriminant;
    ro = x3;
    errorFlag = false;
    iterFlag = false;

    while((iter < m_iter) && (error > max_error_perc))
        x3 = x2 - discriminant;
        ro = x3;
        error = abs(100*(x3 - x2)/x2);
        errors(i) = error;

        if(iter == m_iter)
           iterFlag = true;
        end
        if(error < max_error_perc)
           errorFlag = true;
        end

        iter = iter + 1;
        i = i + 1;
        x0 = x1;
        x1 = x2;
        x2 = x3;

        partialdifferentials0 = (func(x1) - func(x0))/(x1 - x0);
        partialdifferentials1 = (func(x2) - func(x1))/(x2 - x1);
        partial1 = (x2 - x1);
        partial0 = (x1 - x0);
        c = func(x2);
        a = (partialdifferentials1 - partialdifferentials0)/(partial0 + partial1);
        b = a*partial1 + partialdifferentials1;
        discriminant = 2*c/(b + sqrt(b*b - 4*a*c));

    end

    fprintf("The first root of the equation is: %f\n", ro);

    if(errorFlag == true)
        disp("Exceution of the program was stopped since the maximum error criteria as given in the question was met");
    end
    if(iterFlag == true)
        disp("Execution of the program was stopped since the maximum iterations as given in the question were completed");
    end

    figure(1);
    plot(errors);
    title("Plot of Errors vs Iterations");
    xlabel("Iterations");
    ylabel("Relative approximate error");
    grid on;
    figure(2);
    fplot(func);
    title("Plot for the function");
    xlabel("x");
    ylabel("f(x)");
    grid on;

   
% Bairstow Method
elseif(method == '2')
    func = input('What is the polynomial for which we have to find: ','s');
        f = str2sym(func);
        r = input('What is the starting value of r: ');
        s = input('What is the starting value of s: ');
        relativeApproximateError = input('What is the allowed Value of relative error: ');
        max_iterations = input('What is the allowed number of  maximum iteration: ');
        y = matlabFunction(f);
        n = polynomialDegree(f) + 1;
        m = n - 1;
        a = sym2poly(f);
        a = fliplr(a);
        M = zeros;
        N = zeros;
        while n > 0
            for j = 1: max_iterations
                b = zeros;
                c = zeros;
                b(n) = a(n);
                b(n-1) = a(n-1) + r*b(n);
                c(n) = b(n);
                c(n-1) = b(n-1) + r*c(n);

                for i = n-2: -1: 1 
                    b(i) = a(i) + r*b(i+1) + s*b(i+2);
                    c(i) = b(i) + r*c(i+1) + s*c(i+2);
                end

                syms x y 
                eqn1 = c(2)*x + c(3)*y == -b(1);
                eqn2 = c(3)*x + c(4)*y == -b(2);              
                sol = solve([eqn1, eqn2], [x, y]);
                delr = sol.x;
                dels = sol.y;       
                errorPercentr = abs((delr/(r+delr))*100);       
                errorPercents = abs((dels/(r+dels))*100);
                N(j) = errorPercents;
                M(j) = errorPercentr;                
                if errorPercents > relativeApproximateError || errorPercentr > relativeApproximateError 
                    r = r + delr;
                    s = s + dels;
                end 
                if errorPercents < relativeApproximateError && errorPercentr < relativeApproximateError
                    r = r + delr;
                    s = s + dels;
                break;
                end
            end
            m = m - 2;
            if (r^2 + 4*s) < 0
                x1 = r/2;
                x2 = (-(r^2 + 4*s))/2;
                fprintf('Roots of the functions are: %f + %fi, %f - %fi\n', x1, x2, x1, x2);
            else    
                x1 = (r + (r^2 + 4*s)^(1/2))/2;
                x2 = (r - (r^2 + 4*s)^(1/2))/2;    
                fprintf('Roots of the functions are: %f, %f\n', x1, x2);
            end 
            if m > 2
                a = ones;
                for k = n : -1: 3
                   a(k-2) = b(k);
                end  
            end
            n = n - 2;
            if m == 2
                if (b(4)^2 - 4*b(5)*b(3)) < 0
                    x3 = -b(4)/2*b(5);
                    x4 = (-(b(4)^2 - 4*b(5)*b(3)))^(1/2)/2*b(5);
                    fprintf('Roots of the functions are: %f + %fi, %f - %fi\n', x3, x4, x3, x4);
                else
                    x3 = (-b(4) + (b(4)^2 - 4*b(5)*b(3))^(1/2))/2*b(5);
                    x4 = (-b(4) - (b(4)^2 - 4*b(5)*b(3))^(1/2))/2*b(5);
                    fprintf('Roots of the functions are: %f, %f\n', x3, x4);
                end
            break;
            end
            if m == 1
                x1 = -b(3)/b(4);
                fprintf('Roots of the polynomial are: %f\n', x1);
            break;
            end
        end
        figure(1);       
        ezplot(f)
        grid on
        figure(2);
        plot(M)
        title('Error plot of r');
        figure(3);
        plot(N)
        title('Error plot of s');
        
        
    
end