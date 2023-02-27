clear
close all
clc

% program starts and an input is taken from the user for the selection of
% the method 
method = input('Choose the number for the selection of  method of solution :\n Bisection-1,\n False Position-2,\n Fixed Point Method-3,\n Newton-Raphson-4,\n Secant Method-5 \n','s');

% The code for the Bisection Method
if(method == '1')
    clear
    format long
    funcInput = input("Enter the function\n", 's');
    funcInput = strcat('@(x)', funcInput);
    func = str2func(funcInput);
    x1 = input("Enter the  lower limit of the  bracket  for the root  :\n");
    x2 = input("Enter the upper limit of the bracket for the root  :\n");
    disp("Now Enter the stopping criteria");
    m_iter = input("What is the allowed  maximum number of iterations for the algorithm:\n");
    m_error = input("What is the maximum percentage error allowed:\n");

  if x1>=x2
    fprintf('Please enter proper bracket\n');
    return;
  end

    iterations = 0;
    i = 1;
    init = 0;
    mid = (x1 + x2)/2;
    error = abs(100*(mid - init)/mid);
    errorFlag = false;
    iterFlag = false;
    root = mid;

    while((error > m_error) && (iterations < m_iter))
        if(func(x1)*func(x2) < 0)
            mid = (x1 + x2)/2;
            if(func(x1)*func(mid) < 0) 
                x2 = mid;
            else
                x1 = mid;
            end
        end
        error = abs(100*(mid - init)/mid);
        errors(i)=error;
        init = mid;
        root = mid;
        i = i + 1;
        iterations = iterations + 1;
        if(error < m_error)
            errorFlag = true;
        end
        if(iterations == m_iter)
            iterFlag = true;
        end          
    end

    fprintf("The final root of the equation is: %f\n",root);

    if(errorFlag == true)
        disp("Exceution of the program was stopped since the maximum error criteria as given in the question was met");
    end
    if(iterFlag == true)
        disp("Execution of the program was stopped since the maximum iterations as given in the question were completed");
    end

    figure(1);
    plot(errors);
    title("Plot of Errors vs Iterations in Bisection method ");
    xlabel("Iterations");
    ylabel("Relative approximate error");
    grid on;
    figure(2);
    fplot(func);
    title("Plot for the function values");
    xlabel("x");
    ylabel("f(x)");
    grid on;


%% The code for regular falsi method 
elseif(method  == '2')
    clear
    format long
    funcInput = input("Enter the function\n", 's');
    funcInput = strcat('@(x)', funcInput);
    func = str2func(funcInput);
    a0 = input("Enter the lower limit of the bracket for the root :\n");
    b0 = input("Enter the upper limit of the bracket for the root :\n");
  if a0>=b0
    fprintf('Please enter proper bracket\n');
    return;
  end
    disp("Now Enter the stopping criteria");
    m_iter = input("What is the maximum number of iterations allowed for the algorithm:\n");
    m_error = input("What is the maximum percentage error allowed:\n");

    iterations = 0;
    i = 1;
    init = 0;
    cfp = a0 - func(a0)*(b0 - a0)/(func(b0) - func(a0));
    error = abs(100*(cfp - init)/cfp);
    root = cfp;
    errorFlag = false;
    iterFlag = false;

    while((iterations < m_iter) && (error > m_error))
        if(func(a0)*func(b0) < 0)
            cfp = a0 - func(a0)*(b0 - a0)/(func(b0) - func(a0));
            if(func(a0)*func(cfp) < 0)
                b0 = cfp;
            else
                a0 = cfp;
            end
        end

        error = abs(100*(cfp - init)/cfp);
        errors(i) = error;
        init = cfp;
        root = cfp;
        iterations = iterations + 1;
        i = i + 1;
        if(error < m_error)
            errorFlag = true;
        end
        if(iterations == m_iter)
            iterFlag = true;
        end
    end
    
    fprintf("The final root of the equation is: %f\n", root);

    if(errorFlag == true)
        disp("Exceution of the program was stopped since the maximum error criteria as given in the question was met");
    end
    if(iterFlag == true)
        disp("Execution of the program was stopped since the maximum iterations as given in the question were completed");
    end

    figure(1);
    plot(errors);
    title("Plot of Errors vs Iterations in Regular falsi method ");
    xlabel("Iterations");
    ylabel("Relative approximate error");
    grid on;
    figure(2);
    fplot(func);
    title("Plot for the function");
    xlabel("x");
    ylabel("f(x)");
    grid on;

%% Open method
%% Fixed Point Method
elseif(method  == '3')
    clear
    format long
    funcInput = input("Enter the function given in the question\n", 's');
    funcInput = strcat('@(x)', funcInput);
    func = str2func(funcInput);
    phixInput = input("Enter the second function given in the question \n", 's');
    phixInput = strcat('@(x)', phixInput);
    phix = str2func(phixInput);
    a0 = input("Enter the first starting point for initiating the algorithm:\n");
    disp("Now Enter the stopping criteria");
    m_iter = input("What is the maximum number of iterations allowed for the algorithm:\n");
    m_error = input("What is the maximum percentage error allowed:\n");

    iterations = 0;
    i = 1;
    init = 0;
    nextvalue = phix(a0);
    error = abs(100*(nextvalue - a0)/nextvalue);
    root = nextvalue;
    errorFlag = false;
    iterFlag = false;

    while((iterations < m_iter) && (error > m_error))
        nextvalue = phix(a0);
        error = abs(100*(nextvalue - a0)/nextvalue);
        root = nextvalue;

        errors(i) = error;
        a0 = nextvalue;
        iterations = iterations + 1;
        i = i + 1;
        if(error < m_error)
            errorFlag = true;
        end
        if(iterations == m_iter)
            iterFlag = true;
        end
    end
    
    fprintf("The final root of the equation is: %f\n", root);

    if(errorFlag == true)
        disp("Exceution of the program was stopped since the maximum error criteria as given in the question was met");
    end
    if(iterFlag == true)
        disp("Execution of the program was stopped since the maximum iterations as given in the question were completed");
    end

    figure(1);
    plot(errors);
    title("Plot of Errors vs Iterations in the fixed point method");
    xlabel("Iterations");
    ylabel("Relative approximate error");
    grid on;
    figure(2);
    fplot(func);
    title("Plot for the function");
    xlabel("x");
    ylabel("f(x)");
    grid on;

% Newton-Raphson Method
elseif(method  == '4')
    clear 
    format long
    funcInput = input("Enter the function\n", 's');
    derivativeoffunctionInput = input("Enter the first derivative of the function\n", 's');
    funcInput = strcat("@(x)", funcInput);
    derivativeoffunctionInput = strcat("@(x)", derivativeoffunctionInput);
    func = str2func(funcInput);
    funcderivative = str2func(derivativeoffunctionInput);
    xs = input("Enter the starting value for initiating the algorithm:\n");
    disp("Now Enter the stopping criteria");
    m_iter = input("What is the maximum number of iterations allowed for the algorithm:\n");
    m_error = input("What is the maximum percentage error allowed:\n");

    iterations = 0;
    i = 1;
    init = 0;
    nextvalue = xs - func(xs)/funcderivative(xs);
    root = nextvalue;
    error = abs(100*(nextvalue - xs)/nextvalue);
    errorFlag = false;
    iterFlag = false;

    while((iterations < m_iter) && (error > m_error))
       nextvalue = xs - func(xs)/funcderivative(xs);
       error = abs(100*(nextvalue - xs)/nextvalue);

       if(iterations == m_iter)
           iterFlag = true;
       end
       if(error < m_error)
           errorFlag = true;
       end

       errors(i) = error;
       iterations = iterations + 1;
       i = i + 1;
       xs = nextvalue;
       root = nextvalue;

    end

    fprintf("The final root of the equation is: %f\n", root);

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
    ylabel("Approximate Percentage Error");
    grid on;
    figure(2);
    fplot(func);
    title("Plot for the function");
    xlabel("x");
    ylabel("f(x)");
    grid on;

% Secant Method
elseif(method  == '5')
    clear 
    format long
    funcInput = input("Enter the function\n", 's');
    funcInput = strcat('@(x)', funcInput);
    func = str2func(funcInput);
    a0secant  = input("Enter the first starting point for initiating the algorithm:\n");
    b0secant = input("Enter the second starting point for initiating the algorithm:\n");
    disp("Now Enter the stopping criteria");
    m_iter = input("Enter the maximum number of iterations for the algorithm:\n");
    m_error = input("Enter the maximum percentage error allowed:\n");

    iterations = 0;
    i = 1;
    init = 0;
    nextvalue = a0secant - func(a0secant)*(b0secant - a0secant)/(func(b0secant) - func(a0secant));
    root = nextvalue;
    error = 100;
    errorFlag = false;
    iterFlag = false;

    while((iterations < m_iter) && (error > m_error))
       nextvalue = a0secant - func(a0secant)*(b0secant - a0secant)/(func(b0secant) - func(a0secant));
       error = abs(100*(nextvalue - b0secant)/nextvalue);
       root = nextvalue;
       if(iterations == m_iter)
           iterFlag = true;
       end
       if(error < m_error)
           errorFlag = true;
       end

       errors(i) = error;
       iterations = iterations + 1;
       i = i + 1;
       a0secant = b0secant;
       b0secant = nextvalue;

    end

    fprintf("The final root of the equation is: %f\n", root);

    if(errorFlag == true)
        disp("Exceution of the program was stopped since the maximum error criteria as given in the question was met");
    end
    if(iterFlag == true)
        disp("Execution of the program was stopped since the maximum iterations as given in the question were completed");
    end

    figure(1);
    plot(errors);
    title("Plot of Errors vs Iterations for secant method ");
    xlabel("Iterations");
    ylabel("Relative approximate error");
    grid on;
    figure(2);
    fplot(func);
    title("Plot for the function");
    xlabel("x");
    ylabel("f(x)");
    grid on;

end

