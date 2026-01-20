function [t, zsol, dzdt_eval] = odeExplicitSolvers(odefun, tspan, tau, incond, options)

arguments
    odefun function_handle
    tspan (1,2) double
    tau (1,1) double
    incond 
    options.Method (1,1) string {mustBeMember(options.Method, ["RK3", "RK4", "RKB5", "RKN5", "RKB6", "RKB7", "RKCV8"])} = "RK4"
    options.SafeRegime (1,1) logical = false 
end

switch options.Method
    case "RK3"
        c_vector = [0 1/2 1] .';

        A_matrix = [0 0 0;            
                    1/2 0 0;            
                    -1 2 0];
        
        b_vector = [1/6 2/3 1/6] .';
    case "RK4"
        c_vector = [0 1/2 1/2 1] .'; 

        A_matrix = [0 0 0 0; 
                    1/2 0 0 0; 
                    0 1/2 0 0; 
                    0 0 1 0];
        
        b_vector = [1/6 1/3 1/3 1/6] .';
        
    case "RKB5"
        c_vector = [0 1/4 1/4 1/2 3/4 1] .';

        A_matrix = [0 0 0 0 0 0;
                    1/4 0 0 0 0 0;
                    1/8 1/8 0 0 0 0;
                    0 0 1/2 0 0 0;
                    3/16 -3/8 3/8 9/16 0 0;
                    -3/7 8/7 6/7 -12/7 8/7 0];
        
        b_vector = [7/90 0 32/90 12/90 32/90 7/90] .';
        
    case "RKN5"
        c_vector = [0 1/3 2/5 1 2/3 4/5] .';

        A_matrix = [0 0 0 0 0 0;
                    1/3 0 0 0 0 0;
                    4/25 6/25 0 0 0 0;
                    1/4 -3 15/4 0 0 0;
                    2/27 10/9 -50/81 8/81 0 0;
                    2/25 12/25 2/15 8/75 0 0];
        
        b_vector = [23/192 0 125/192 0 -27/64 125/192] .';

    case "RKB6"
        c_vector = [0 1/3 2/3 1/3 5/6 1/6 1] .';

        A_matrix = [0 0 0 0 0 0 0;
                    1/3 0 0 0 0 0 0;
                    0 2/3  0 0 0 0 0;
                    1/12 1/3 -1/12 0 0 0 0;
                    25/48 -55/24 35/48 15/8  0 0 0;
                    3/20 -11/24 -1/8 1/2 1/10 0 0;
                    -261/260 33/13 43/156 -118/39 32/195 80/39 0];
  
        b_vector = [13/200 0 11/40 11/40 4/25 4/25 13/200] .';

    case "RKB7"
        c_vector = [0 1/6 1/3 1/2 2/11 2/3 6/7 0 1] .';

        A_matrix = [0 0 0 0 0 0 0 0 0;
                    1/6 0 0 0 0 0 0 0 0;
                    0 1/3 0 0 0 0 0 0 0;
                    1/8 0 3/8 0 0 0 0 0 0;
                    148/1331 0 150/1331 -56/1331 0 0 0 0 0;
                    -404/243 0 -170/27 4024/1701 10648/1701 0 0 0 0;
                    2466/2401 0 1242/343 -19176/16807 -51909/16807 1053/2401 0 0 0;
                    5/154 0 0 96/539 -1815/20384 -405/2464 49/1144 0 0;
                    -113/32 0 -195/22 32/7 29403/3584 -729/512 1029/1408 21/16 0];
        
        b_vector = [0 0 0 32/105 1771561/6289920 243/2560 16807/74880 77/1440 11/270] .';

    case "RKCV8"
        c_vector = [0 1/2 1/2 (7+sqrt(21))/14 (7+sqrt(21))/14 1/2 (7-sqrt(21))/14 (7-sqrt(21))/14 1/2 (7+sqrt(21))/14 1] .';

        A_matrix = [0 0 0 0 0 0 0 0 0 0 0;
                    1/2 0 0 0 0 0 0 0 0 0 0;
                    1/4 1/4 0 0 0 0 0 0 0 0 0;
                    1/7 (-7-3*sqrt(21))/98 (21+5*sqrt(21))/49 0 0 0 0 0 0 0 0;
                    (11+sqrt(21))/84 0 (18+4*sqrt(21))/63 (21-sqrt(21))/252 0 0 0 0 0 0 0;
                    (5+sqrt(21))/48 0 (9+sqrt(21))/36 (-231+14*sqrt(21))/360 (63-7*sqrt(21))/80 0 0 0 0 0 0;
                    (10-sqrt(21))/42 0 (-432+92*sqrt(21))/315 (633-145*sqrt(21))/90 (-504+115*sqrt(21))/70 (63-13*sqrt(21))/35 0 0 0 0 0;
                    1/14 0 0 0 (14-3*sqrt(21))/126 (13-3*sqrt(21))/63 1/9 0 0 0 0;
                    1/32 0 0 0 (91-21*sqrt(21))/576 11/72 (-385-75*sqrt(21))/1152 (63+13*sqrt(21))/128 0 0 0;
                    1/14 0 0 0 1/9 (-733-147*sqrt(21))/2205 (515+111*sqrt(21))/504 (-51-11*sqrt(21))/56 (132+28*sqrt(21))/245 0 0;
                    0 0 0 0 (-42+7*sqrt(21))/18 (-18+28*sqrt(21))/45 (-273-53*sqrt(21))/72 (301+53*sqrt(21))/72 (28-28*sqrt(21))/45 (49-7*sqrt(21))/18 0];
        
        b_vector = [1/20 0 0 0 0 0 0 49/180 16/45 49/180 1/20] .';
end

c = c_vector(:);
A = A_matrix;
b = b_vector(:);
z_0 = incond(:);
dzdt = odefun;

s_stages = numel(c);
m = numel(z_0);

t = (tspan(1) : tau : tspan(2)) .';
if t(end) ~= tspan(2)
    if (tau > 0 && t(end) < tspan(2)) || (tau < 0 && t(end) > tspan(2))
        t(end+1) = tspan(2);
    elseif abs(t(end) - tspan(2)) < tol
        t(end) = tspan(2);
    end
end

N = numel(t);
zsol = zeros(m, N);
K = zeros(m, s_stages);
dzdt_eval = zeros(N, m);

zsol(:, 1) = z_0;

if options.SafeRegime
    for n = 1:N - 1
        tau_n = t(n+1) - t(n);
        K(:, 1) = dzdt(t(n), zsol(:, n)); 
        
        dzdt_eval(n, :) = K(:, 1) .';
    
        for i = 2:s_stages
            K(:, i) = dzdt( ...
                t(n) + tau_n * c(i), ...
                zsol(:, n) + tau_n * K(:, 1:i-1) * A(i, 1:i-1).' ...
            );
        end
    
        zsol(:, n+1) = zsol(:, n) + tau_n * K * b;
    
        if any(~isfinite(zsol(:, n+1)))
            warning('Integration stopped: solution contains NaN or Inf at step %d (t = %g).', n+1, t(n+1));
            t = t(1:n+1);
            zsol = zsol(:, 1:n+1);
            dzdt_eval = dzdt_eval(1:n+1, :);
            break;
        end
    end
else
    for n = 1:N - 1
        tau_n = t(n+1) - t(n);
        K(:, 1) = dzdt(t(n), zsol(:, n)); 
        
        dzdt_eval(n, :) = K(:, 1) .';
    
        for i = 2:s_stages
            K(:, i) = dzdt( ...
                t(n) + tau_n * c(i), ...
                zsol(:, n) + tau_n * K(:, 1:i-1) * A(i, 1:i-1).' ...
            );
        end
    
        zsol(:, n+1) = zsol(:, n) + tau_n * K * b;
    end    
end

dzdt_eval(end, :) = dzdt(t(end), zsol(:, end)) .';

zsol = zsol .';

if (nargout < 3)
    dzdt_eval = [];
end

end