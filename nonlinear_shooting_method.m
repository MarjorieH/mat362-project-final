% endpoints a & b
% boundary conditions alpha and beta
% number of subintervals N (>=2)
% tolerance TOL
% maximum iterations M
function [xv, W] = nonlinear_shooting_method(f, fy, fyp, a, b, alpha, beta, N, TK, TOL, M)

    h = (b - a) / N; % step size
    k = 1; % iterator
    %TK = (beta - alpha) / (b - a); % slope
    xv = 1:N;
    
    while (k <= M)
        
        W = zeros([2, N]);
        K = zeros([4, 2]);
        Kp = zeros([4, 2]);
                
        W(1,1) = alpha;
        W(2,1) = TK;
        u1 = 0;
        u2 = 1;
        
        for i = 2:N
                        
            x = a + (i - 1) * h;
            
            K(1,1) = h * W(2,i-1);
            K(1,2) = h * f(x, W(1,i-1), W(2,i-1));
            K(2,1) = h * (W(2,i-1) + K(1,2) / 2);
            K(2,2) = h * f(x + h / 2, W(1,i-1) + K(1,1) / 2, W(2,i-1) + K(1,2) / 2);
            K(3,1) = h * (W(2,i-1) + K(2,2) / 2);
            K(3,2) = h * f(x + h / 2, W(1,i-1) + K(2,1) / 2, W(2,i-1) + K(2,2) / 2);
            K(4,1) = h * (W(2,i-1) + K(3,2));
            K(4,2) = h * f(x + h, W(1,i-1) + K(3,1), W(2,i-1) + K(3,2));
            
            W(1,i) = W(1,i-1) + (K(1,1) + 2 * K(2,1) + 2 * K(3,1) + K(4,1)) / 6;
            W(2,i) = W(2,i-1) + (K(1,2) + 2 * K(2,2) + 2 * K(3,2) + K(4,2)) / 6;
            
            Kp(1,1) = h * u2;
            Kp(1,2) = h * (fy(x, W(1,i-1), W(2,i-1)) * u1 + fyp(x, W(1,i-1), W(2,i-1)) * u2);
            Kp(2,1) = h * (u2 + Kp(1,2) / 2);
            Kp(2,2) = h * (fy(x + h / 2, W(1,i-1), W(2,i-1)) * (u1 + Kp(1,1) / 2) + fyp(x + h / 2, W(1,i-1), W(2,i-1)) * (u2 + Kp(1,2) / 2));
            Kp(3,1) = h * (u2 + Kp(2,2) / 2);
            Kp(3,2) = h * (fy(x + h / 2, W(1,i-1), W(2,i-1)) * (u1 + Kp(2,1) / 2) + fyp(x + h / 2, W(1,i-1), W(2,i-1)) * (u2 + Kp(2,2) / 2));
            Kp(4,1) = h * (u2 + Kp(3,2));
            Kp(4,2) = h * (fy(x + h / 2, W(1,i-1), W(2,i-1)) * (u1 + Kp(3,1) / 2) + fyp(x + h / 2, W(1,i-1), W(2,i-1)) * (u2 + Kp(3,2) / 2));
            
            u1 = u1 + (Kp(1,1) + 2 * Kp(2,1) + 2 * Kp(3,1) + Kp(4,1)) / 6;
            u2 = u2 + (Kp(1,2) + 2 * Kp(2,2) + 2 * Kp(3,2) + Kp(4,2)) / 6;
            
        end
        
        if abs(W(1,N) - beta) <= TOL
            
            for i = 1:N
                xv(i) = a + i * h;
            end
            return;
        end
        
        TK = TK - (W(1,N) - beta) / u1;
        k = k + 1;
    end
    
    fprintf('Maximum number of iterations exceeded.');

end

