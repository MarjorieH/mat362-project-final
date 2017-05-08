a = 0;
b = 1;
N = 100;

c = 0;
d = 1;
eqn1 = @(y, yp) yp;
eqn2 = @(y, yp) -y.^3;
t = linspace(a, b, N+1)';

diff = 1;
while diff > 0
    [y1, y2] = rk_mod(eqn1, eqn2, t, c, d);
    plot(t, y1);
    hold on;
    d = d + 0.1; % increase d slightly
    diff = y1(N+1);
    % print in table format for LaTex
    fprintf('%g & %f \\\\ \\hline\n', d, y1(N+1));
end
