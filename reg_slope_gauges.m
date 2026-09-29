clc; clear;

% Gauge
fdir = 'C:\Users\User\Desktop\FUNWAVE-TVD-Version_3.6\FUNWAVE-TVD-Version_3.6\simple_cases\surface_wave_1d\reg_slope\output\';
x_locations = [400, 600, 1000];
x_wm = 250;
a = 0.5;
h = 10;
T = 8;
g = 9.81;
w=2*pi/T;

% Newton-Raphson methon

f = @(k) g*k*tanh(k*h)-w^2;
df = @(k) g*tanh(k*h)+g*k*h*(sech(k*h))^2;

x1 = 0.0001;  % 下限
x2 = 1.0;     % 上限

k = w^2*h/g;

for i = 1 : 5000;
    fk = f(k);
    dfk = df(k);

    fprintf('%4d     %.8f     %.10f\n', i, k, fk);

    
    if abs(fk) < 10^-5;
        break;
    end

    k_n = k - fk / dfk;
    
    if (k_n > x1 && k_n < x2);
        k = k_n;
    else 
        fprintf('error\n');
    end
end

L = 2*pi/k;
fprintf('週波數 k = %.6f (rad/m)\n', k);
fprintf('波長   L = %.6f (m)\n', L);

figure(3); clf;

for n = 1:3
    fname = sprintf('sta_000%d', n);
    data = importdata([fdir fname]);
    [t_sim, idx] = sort(data(:, 1));
    eta_sim = data(idx, 2);
    
    % 理論值
    dist = x_locations(n) - x_wm;
    eta_theory = a * cos(k * dist - w * t_sim + 4.8);
    
    % 到達時間
    Cg = (w/k) * 0.5 * (1 + (2*k*h)/sinh(2*k*h));
    eta_theory(t_sim < dist/Cg) = 0;
    
    % 繪圖
    subplot(3, 1, n);
    plot(t_sim, eta_sim, 'b', 'LineWidth', 1.2); hold on;
    plot(t_sim, eta_theory, 'r--', 'LineWidth', 1.0);
    
    grid on;
    ylabel('$\eta$ (m)', 'Interpreter', 'latex');
    title(['Gauge at x = ', num2str(x_locations(n)), ' m']);
    xlim([0 200]); 
    ylim([-1 1]);
    
    if n == 1
        legend('FUNWAVE','Linear Theory','Location','northeast');
    end
end
xlabel('Time (s)');