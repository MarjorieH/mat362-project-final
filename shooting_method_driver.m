A = 0;
B = 50;
N = 100;

s = linspace(A, B, N+1);

D = 1:N+1;

last = 0; % saves the index of the last eigenvalue
figure();

for i=1:N+1

    a = 0;
    b = 1;
    n = 60;

    eqn1 = @(y, yp) yp;
    eqn2 = @(y, yp) -y.^3 - s(i) .* y;
    t = linspace(a, b, n+1)';
    
    for d=0:0.1:200
        
        % solve the IVP
        [y1, y2] = rk_mod(eqn1, eqn2, t, 0, d);
        
        % if the result solves the BVP, save the solution
        if y1(n) * y1(n+1) < 0
            D(i) = d;
            break;
        end
    end
    
    % if we hit an eigenvalue, plot the branch for the eigenfunction
    if abs(D(i)) <= 0.1
        plot(s(last+1:i), D(last+1:i));
        hold on;
        last = i;
    end
end
% plot the remaining branch as well, if needed
plot(s(last+1:end),D(last+1:end));


