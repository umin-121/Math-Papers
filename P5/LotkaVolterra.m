

clear; clc; close all;


%  x' = x(a1 - b11 x - b12 y),  y' = y(-a2 + b21 x - b22 y)
b = [1 1; 1 1];                 % [b11 b12; b21 b22]

% Case A: a1/b11 > a2/b21
aA = [4; 1];
[xA, yA] = equil2(aA, b, 0, 0);
fprintf('--- Section 2, case A: a = [4 1] ---\n');
fprintf('Interior equilibrium (x,y) = (%.4f, %.4f)   [paper: (2.5, 1.5)]\n', xA, yA);
report2(aA, b, 0, 0, xA, yA);
plotPhase(aA, b, 0, 0, 'System (1), a_1/b_{11} > a_2/b_{21}');

% Case B : a1/b11 < a2/b21
aB = [1; 4];
[xB, yB] = equil2(aB, b, 0, 0);
fprintf('\n--- Section 2, case B: a = [1 4] ---\n');
fprintf('Nullcline intersection = (%.4f, %.4f): y < 0, outside first quadrant\n', xB, yB);
fprintf('Stable equilibrium is the prey-only point (%.4f, 0): predators die out\n', aB(1)/b(1,1));
report2(aB, b, 0, 0, aB(1)/b(1,1), 0);
plotPhase(aB, b, 0, 0, 'System (1), a_1/b_{11} < a_2/b_{21}');


%  x' = x(a1 - b11 x - b12 y) - alpha x
%  y' = y(-a2 + b21 x - b22 y) - beta y

alpha = 0.3;  beta = 0.2;
[xH, yH] = equil2(aA, b, alpha, beta);
fprintf('\n--- Section 3: alpha = %.1f, beta = %.1f ---\n', alpha, beta);
fprintf('Interior equilibrium (x,y) = (%.4f, %.4f)   [paper: (2.45, 1.25)]\n', xH, yH);
report2(aA, b, alpha, beta, xH, yH);

mWar   = yA / xA;    % wartime slope  (no fishing)
mPeace = yH / xH;    % peacetime slope (fishing)
fprintf('Wartime slope   m_a = %.6f  (paper: 3/5   = %.6f)\n', mWar,   3/5);
fprintf('Peacetime slope m_b = %.6f  (paper: 25/49 = %.6f)\n', mPeace, 25/49);
fprintf('Harvesting shifts the equilibrium: m_b < m_a is %s\n', mat2str(mPeace < mWar));

plotPhase(aA, b, alpha, beta, 'System (2) with harvesting');


%  x = scallops, y = rays, z = sharks
%  x' = x(a1 - b11 x - b12 y) - alpha x
%  y' = y(a2 - b22 y + b21 x - b23 z) - beta y
%  z' = z(-a3 - b33 z + b32 y) - gamma z

a3s = [4; 2; 0.5];                       % [a1; a2; a3]
B3  = [1 1 0;                            % [b11 b12  0 ]
       1 1 1;                            % [b21 b22 b23]
       0 1 1];                           % [ 0  b32 b33]
alpha3 = 0.3;  beta3 = 0.2;  gamma0 = 0;

[eq0, M0] = equil3(a3s, B3, alpha3, beta3, gamma0);
fprintf('\n--- Section 4: three species, gamma = %.1f ---\n', gamma0);
fprintf('Interior equilibrium (x,y,z) = (%.4f, %.4f, %.4f)\n', eq0);
fprintf('Eigenvalue real parts of Jacobian: %s\n', mat2str(real(eig(jac3(eq0,a3s,B3,alpha3,beta3,gamma0))).', 4));
fprintf('det(M) = %.4f\n', det(M0));

% Sensitivity of the equilibrium to shark harvesting )
h = 1e-5;
eqp = equil3(a3s, B3, alpha3, beta3, gamma0 + h);
eqm = equil3(a3s, B3, alpha3, beta3, gamma0 - h);
d = (eqp - eqm) / (2*h);
fprintf('dx/dgamma = %+.4f (scallops)\n', d(1));
fprintf('dy/dgamma = %+.4f (rays)\n',     d(2));
fprintf('dz/dgamma = %+.4f (sharks)\n',   d(3));



function du = sys2(~, u, a, b, alpha, beta)
    x = u(1); y = u(2);
    du = [ x*(a(1) - b(1,1)*x - b(1,2)*y) - alpha*x;
           y*(-a(2) + b(2,1)*x - b(2,2)*y) - beta*y ];
end

function [x, y] = equil2(a, b, alpha, beta)
    %   b11 x + b12 y =  a1 - alpha
    %  -b21 x + b22 y = -(a2 + beta)
    M = [ b(1,1)  b(1,2);
         -b(2,1)  b(2,2) ];
    r = [ a(1) - alpha;
         -(a(2) + beta) ];
    s = M \ r;
    x = s(1); y = s(2);
end

function report2(a, b, alpha, beta, x, y)
    J = [ a(1)-alpha-2*b(1,1)*x-b(1,2)*y,      -b(1,2)*x;
          b(2,1)*y,   -(a(2)+beta)+b(2,1)*x-2*b(2,2)*y ];
    fprintf('Jacobian eigenvalues at (%.3f, %.3f): %s\n', x, y, mat2str(eig(J).', 4));
end

function plotPhase(a, b, alpha, beta, ttl)
    figure; hold on; box on
    [X0, Y0] = meshgrid(0.25:0.75:4.75, 0.25:0.75:4.75);
    for k = 1:numel(X0)
        [~, U] = ode45(@(t,u) sys2(t,u,a,b,alpha,beta), [0 60], [X0(k); Y0(k)]);
        plot(U(:,1), U(:,2), 'Color', [0.65 0.75 0.95]);
    end
    xx = linspace(0, 5, 200);
    h1 = plot(xx, (a(1)-alpha - b(1,1)*xx)/b(1,2), 'g-', 'LineWidth', 1.8);
    h2 = plot(xx, (-(a(2)+beta) + b(2,1)*xx)/b(2,2), 'r-', 'LineWidth', 1.8);
    [xe, ye] = equil2(a, b, alpha, beta);
    h3 = plot(0, 0, 'ko', 'MarkerFaceColor', 'k');
    plot((a(1)-alpha)/b(1,1), 0, 'ko', 'MarkerFaceColor', 'k');
    if ye > 0
        plot(xe, ye, 'ko', 'MarkerFaceColor', 'k');
    end
    axis([0 5 0 5]); xlabel('x (prey)'); ylabel('y (predator)');
    title(ttl);
    legend([h1 h2 h3], 'prey nullcline', 'predator nullcline', 'equilibria', ...
           'Location', 'northeast');
end

function [e, M] = equil3(a, B, alpha, beta, gamma)
    % Interior equilibrium: 
    M = [ B(1,1)   B(1,2)   0;
         -B(2,1)   B(2,2)   B(2,3);
          0       -B(3,2)   B(3,3) ];
    r = [ a(1) - alpha;
          a(2) - beta;
         -(a(3) + gamma) ];
    e = (M \ r).';
end

function J = jac3(e, a, B, alpha, beta, gamma)
    x = e(1); y = e(2); z = e(3);
    J = [ a(1)-alpha-2*B(1,1)*x-B(1,2)*y,  -B(1,2)*x,  0;
          B(2,1)*y,  a(2)-beta-2*B(2,2)*y+B(2,1)*x-B(2,3)*z,  -B(2,3)*y;
          0,  B(3,2)*z,  -a(3)-gamma-2*B(3,3)*z+B(3,2)*y ];
end