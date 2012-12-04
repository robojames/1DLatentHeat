clear all; clc;
M1 = xlsread('Data_new_baseline.xls');
M2 = xlsread('Data_new_sameLengthyplus.xls');
M3 = xlsread('Data_05micron.xls');
M4 = xlsread('Data_0p1micron.xls');

t = M1(:,1);
T1 = M1(:,2) - 273;
T2 = M2(:,2) - 273;
T3 = M3(:,2) - 273;
T4 = M4(:,2) - 273;

figure(1); hold on; grid on;
xlabel(' Time, s '); ylabel(' Temperature, C ');
plot(t, T1, 'k*-'); plot(t, T2, '^-'); plot(t, T3, 'k^-'); plot(t, T4);
legend('Baseline', 'Nickel-Chromium - 25e-6 m L_c', 'Nickel-Chromium - 5e-6 m L_c', 'Nickel-Chromium - 0.1e-6 m L_c');

figure(2); hold on; grid on;
xlabel(' Time, s'); ylabel(' Temperature - Baseline, C ');
plot (t, abs(T2-T1), '^-'); plot(t,abs(T3-T1), 'k^-'); plot(t,abs(T4-T1),'r^-');

legend('25e-6m L_c', '5e-6m L_c','0.1e-6m L_c');