% reset
clear variables;
clc;
close all;
 
% khai bao hang so
G = 59292.56; % dv: N
Gtck = 3765.59; % dv: N
Ltck = 5.064; % dv: m
OA = Ltck/2; % dv: m
OM = 0.844; % dv: m
ON = OA-OM; % dv: m
 
% khai bao buoc chay vong lap
CN_start = 0.8;
CN_end = 3.2;
CN_num_steps = 500; % so lan chia phan tu huu han
% tao Array cho CN
CN_step_size = (CN_end - CN_start) / (CN_num_steps - 1);
CN_values = linspace(CN_start, CN_end, CN_num_steps);
 
result_table = zeros(CN_num_steps, 9);
 
for i = 1:CN_num_steps
    % tinh to hop luc Ay+Cy
    AyCy = G;
    CN = CN_values(i);
    
    % tinh goc luong giac
    alpha = asin(CN / Ltck);
    teta1 = 2 * alpha;
    MN = sqrt(OM^2 + ON^2 - 2 * OM * ON * cos(teta1));
    cos_teta2 = (ON^2 + MN^2 - OM^2) / (2 * ON * MN);
    teta2 = acos(cos_teta2);
    beta = alpha + teta2;
    
    AC = sqrt(Ltck^2 - CN^2);
    % tinh force
    NQ = OM;
    MP = OA+OM;
    force_tu = AyCy*Ltck*cos(alpha) + Gtck*(Ltck/2)*cos(alpha);
    force_mau = sin(alpha)*cos(beta)*(NQ+MP-Ltck) + sin(beta)*cos(alpha)*(MP-NQ);
    force = force_tu/force_mau;
    
    % tinh cac phan luc lien ket
    Mx = force*cos(beta);
    My = force*sin(beta);
    CE = 5.0004096 - AC;
    Cy = 2.25*G/(5-CE);
    Ay = AyCy - Cy;
    Ox = -Mx;
    Oy = -Cy*Ltck+My*MP+Mx*MP*tan(alpha)-(Gtck/2)*(Ltck/2)+Ox*(Ltck/2)*tan(alpha);
    Qy = Ay-Oy+My;
    Py = Cy-My+Oy;
    force_A = Ay;
    force_C = Cy;
    force_O = sqrt(Ox^2 + Oy^2);
    force_Q = Qy;
    force_P = Py;
    
    result_table(i, :) = [rad2deg(alpha), force, force_A, force_C, force_O, force_Q, force_P, Ox, Oy];
end
 
% tra ve ket qua
disp('Result Table:');
disp(result_table);

% ve bieu do alpha, force, force_A
figure;
ax = plotyy(result_table(:, 1), result_table(:, 3), result_table(:, 1), result_table(:, 2));
% ten bieu do
ylabel(ax(1), 'Ay (N)');
ylabel(ax(2), 'Fxl (N)');
xlabel('Alpha (degrees)');
% chia nho phan tu, so lan chia: 10
set(ax(1), 'YTick', linspace(min(result_table(:, 3)), max(result_table(:, 3)), 10));
set(ax(2), 'YTick', linspace(min(result_table(:, 2)), max(result_table(:, 2)), 10));
set(ax, 'XTick', linspace(min(result_table(:, 1)), max(result_table(:, 1)), 10));
grid on;

% ve bieu do alpha, force, force_C
figure;
ax = plotyy(result_table(:, 1), result_table(:, 4), result_table(:, 1), result_table(:, 2));
ylabel(ax(1), 'Cy (N)');
ylabel(ax(2), 'Fxl (N)');
xlabel('Alpha (degrees)');
set(ax(1), 'YTick', linspace(min(result_table(:, 4)), max(result_table(:, 4)), 10));
set(ax(2), 'YTick', linspace(min(result_table(:, 2)), max(result_table(:, 2)), 10));
set(ax, 'XTick', linspace(min(result_table(:, 1)), max(result_table(:, 1)), 10));
grid on;

% ve bieu do alpha, force, force_O
figure;
ax = plotyy(result_table(:, 1), result_table(:, 5), result_table(:, 1), result_table(:, 2));
ylabel(ax(1), 'sqrt(Ox^2 + Oy^2) (N)');
ylabel(ax(2), 'Fxl (N)');
xlabel('Alpha (degrees)');
set(ax(1), 'YTick', linspace(min(result_table(:, 5)), max(result_table(:, 5)), 10));
set(ax(2), 'YTick', linspace(min(result_table(:, 2)), max(result_table(:, 2)), 10));
set(ax, 'XTick', linspace(min(result_table(:, 1)), max(result_table(:, 1)), 10));
grid on;

% ve bieu do alpha, force, force_Q
figure;
ax = plotyy(result_table(:, 1), result_table(:, 6), result_table(:, 1), result_table(:, 2));
ylabel(ax(1), 'Qy (N)');
ylabel(ax(2), 'Fxl (N)');
xlabel('Alpha (degrees)');
set(ax(1), 'YTick', linspace(min(result_table(:, 6)), max(result_table(:, 6)), 10));
set(ax(2), 'YTick', linspace(min(result_table(:, 2)), max(result_table(:, 2)), 10));
set(ax, 'XTick', linspace(min(result_table(:, 1)), max(result_table(:, 1)), 10));
grid on;

% ve bieu do alpha, force, force_P
figure;
ax = plotyy(result_table(:, 1), result_table(:, 7), result_table(:, 1), result_table(:, 2));
ylabel(ax(1), 'Py (N)');
ylabel(ax(2), 'Fxl (N)');
xlabel('Alpha (degrees)');
set(ax(1), 'YTick', linspace(min(result_table(:, 7)), max(result_table(:, 7)), 10));
set(ax(2), 'YTick', linspace(min(result_table(:, 2)), max(result_table(:, 2)), 10));
set(ax, 'XTick', linspace(min(result_table(:, 1)), max(result_table(:, 1)), 10));
grid on;