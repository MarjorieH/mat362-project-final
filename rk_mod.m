% rk_mod solves a system of differential equations
% @param f1 - The first differential eqtn
% @param f2 - the second differential eqtn
% @param x - A column vector of step values
% @param c - The initial value for f1
% @param d - The initial value for f2
% Returns - 2 column vectors, y1 and y2, that have the solutions to the
% eqtns f1 and f2 over the interval represented by x
function [y1,y2] = rk_mod(f1, f2, x, c, d)
    % x is a vector of the steps, so size of x will be # of steps to run
    [m,~] = size(x);
    % create vectors to hold the approximated answers for each eqtn
    y1 = ones(m,1);
    y2 = ones(m,1);
    
    % create an intermediate vector to be used in the loop (based on the
    % algorithm in the book)
    w = ones(m,1);
    
    % set the initial values of the only two elements of w we will need -- 
    % we will be overwriting these as we go
    w(1) = c;
    w(2) = d;
    
    % set the initial values of the returned vectors
    y1(1) = c;
    y2(1) = d;
    
    % determine our dx or delta t or h...I like to call it h :)
    h = x(2) - x(1);
    
    % Run from i = 1 to i = m - 1
    for i = 1: m-1
        % do every RK step twice: once for each equation
        k1_1 = h*f1(w(1), w(2));
        k1_2 = h*f2(w(1), w(2));
        
        k2_1 = h*f1(w(1) + k1_1/2, w(2) + k1_2/2);
        k2_2 = h*f2(w(1) + k1_1/2, w(2) + k1_2/2);
        
        k3_1 = h*f1(w(1) + k2_1/2, w(2) + k2_2/2);
        k3_2 = h*f2(w(1) + k2_1/2, w(2) + k2_2/2);
        
        k4_1 = h*f1(w(1) + k3_1, w(2) + k3_2);
        k4_2 = h*f2(w(1) + k3_1, w(2) + k3_2);
        
        w(1) = w(1) + (k1_1 + 2*k2_1 + 2*k3_1 + k4_1)/6;
        w(2) = w(2) + (k1_2 + 2*k2_2 + 2*k3_2 + k4_2)/6;
        
        % store the current answer
        y1(i+1) = w(1);
        y2(i+1) = w(2);
    end
end