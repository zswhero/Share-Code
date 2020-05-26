% code for the one dimension beam propagation
clc 
clear all
% close all

%% Parameters table

w_in = 2; % um, IWs gaussian beam width
w_p = 100; % um, IWs periods
lambda = 1; % um, wavelength 
% zt = 2.*w_p.^2/lambda; % talbot cavity length
zt = 1e3;
k = 2*pi/lambda; % wave vector


% parameters
Ein_test = @(x) exp(-(x./w_in).^2);
Power_Ein = integral(@(x)(Ein_test(x).^2), -Inf, +Inf);


%% Step 1 Input Waveguide (IW) far field 
Ein = @(x, a) ((1/sqrt(Power_Ein)).*exp(-((x - a)./w_in).^2));
E_slit = @(x, a) 1.*(x>=- a./2 & x<= a./2);


x = linspace(-1e2, 1e2, 1e3);
% Example: a gaussian function
Ein_value = Ein(x, 0);
% % Example: a slit function
% Ein_value = E_slit(x, 10);
% % Example: a higher mode function (two gaussian)
% Ein_value = Ein(x, -2) + Ein(x, 2);

x_ = linspace(-3e3, 3e3, 1e5);
k = 2*pi/lambda; 
for i_x_ = 1:length(x_)
    z = zt;
    r0 = sqrt((x_(i_x_) - x).^2 + zt.^2);
    E_core = Ein_value.*exp(1j.*k.*sqrt((x_(i_x_) - x).^2 + zt.^2))./(1j.*lambda.*r0);
    E_f(i_x_) = trapz(x, E_core);
end

% test E_f vs x
figure
power_E_f = trapz(x_, abs(E_f).^2);
plot(x_, (1/sqrt(power_E_f)).*abs(E_f))
xlabel('x_ (um)')
ylabel('|E_f_value|')
% test E_f vs theta
figure
power_E_f = trapz(x_, abs(E_f).^2);
plot(rad2deg( atan(x_./zt)) , (1/sqrt(power_E_f)).*abs(E_f))
xlabel('theta (degree)')
ylabel('|E_f_value|')