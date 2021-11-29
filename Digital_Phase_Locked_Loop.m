% Digital Phase Locked Loop

%% Parameters
Fi = 1e6;               % Input Frequency to the DPLL
Fs = 20e6;              % Sampling Frequency
SpC = Fs/Fi;            % Samples per cycle of the input frequency
Omega_0 = 2*pi/SpC;     % Digital frequency (Radian/sample)
N = 10000;              % Number of samples
eta  = 0.5;             % Loop damping factor
Bn_times_T = 0.001;     % Bn (noise bandwidth) multiplied by the sampling period 
Kp = 1;                 % Phase detector constant
K0 = 1;                 % DDS constant

%% Calculate loop constants

% Theta_n (a definition in terms of T and Bn)
Theta_n = (Bn_times_T)/(eta + (1/(4*eta)));

% Constants obtained by analogy to the continuous time transfer function:
Kp_K0_K1 = (4 * eta * Theta_n) / (1 + 2*eta*Theta_n + Theta_n^2);
Kp_K0_K2 = (4 * Theta_n^2) / (1 + 2*eta*Theta_n + Theta_n^2);

% K1 and K2 (PI constants):
K1 = Kp_K0_K1/(Kp*K0);
K2 = Kp_K0_K2/(Kp*K0);

%% DPLL in Frequency Domain

% Find coeffs of closed-loop transfer function of DDS_phase/input_phase
% Hd(z) = (b1z^-1 + b2z^-2)/(1 + a1Z^-1 + a2z^-2)
b1= Kp*K0*(K1+K2);
b2= -Kp*K0*K1;
a1_a2_constants
b= [0 b1 b2];   % numerator coeffs
a= [1 a1 a2];   % denominator coeffs
% step response
n= 1:N;
t= n/Fs;
x= pi*ones(1,N);        % step function of pi height
%x = pi*(0:N-1)*1/N;    % ramp function
y= filter(b,a,x);       % step or ramp response
pe = x - y;             % phase error response

figure(1)
plot(t*1e3,pe),grid
ylabel('\Delta\theta')
xlabel('ms')
title('Linear Loop Phase Error')

% Closed-loop frequency response
figure(2)
u = 0:.1:.9;
f= 10* 10 .^u;          % log-scale frequencies
f = [f 10*f 100*f 1000*f];
z = exp(1i*2*pi*f/Fs);  % complex frequency z
Hd= (b1*z.^-1 + b2*z.^-2)./(1 + a1*z.^-1 + a2*z.^-2); % closed-loop response
Hd_dB= 20*log10(abs(Hd));
semilogx(f,Hd_dB),grid
xlabel('Hz'),ylabel('Hd(z) dB')

%% DPLL in Time Domain

% Preallocate
Loop_input = zeros(N, 1);
Theta_input = zeros(N, 1);
Theta_output = zeros(N, 1);
Out_cos   = zeros(N, 1);
Out_sin   = zeros(N, 1);
ee   = zeros(N, 1);
vv   = zeros(N, 1);

% Initialize

DDS       = 0;
V21       = 0;
W1        = Omega_0;

for n=1:N
 
      Theta_in = Omega_0*n + pi;
      Theta_input(n) = Theta_in;
      Input = 2*cos(Theta_in);
      Loop_input(n) = Input;
            
      % Evaluate arithmetic expresions in topological order
                
          % Compute dds
             dds = DDS + W1;
             Theta_output(n) = dds;              
                
                
          % Save for plotting:
                
        Out_cos(n)   = 2*cos(dds);
        Out_sin(n)   = sin(dds);
        
        % Phase Error Detector:
                
        e = Loop_input(n)*(-Out_sin(n));
        
        ee(n)   = e;
                                             
        
        %% Loop Filter
    v1   = K1*e;            % Proportional
    v2   = V21 + K2*e;      % Integral
    v    = v1 + v2;         % PI Output
    
    vv(n)   = v;
    
    % Compute DDS input
    w    = Omega_0 + v;      
          
    V21 = v2;
    
    DDS = dds;
    W1 = w;
            
   
end

%% Plots

windowSize = 10; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
ee_filt = filter(b,a,ee);

figure(3)                      % Plot of filtered error (without the second harmonic)
plot(t*1e3,ee_filt),grid
ylabel('e(n)')
xlabel('ms')
title('Filtered Loop Error (Phase Detector Output)')

figure(4)                      % Plot of Loop Phase Error
plot(t*1e3,Theta_input-Theta_output),grid
ylabel('\Delta\theta')
xlabel('ms')
title('Loop Phase Error (\thetain - \thetaout)')

figure(5)                      % Plot of Loop Input and Output over various time ranges
subplot (3,2,1)
plot(t*1e3,Loop_input)
hold on
plot(t*1e3,Out_cos, 'r')
ylabel('Amplitude')
xlabel('ms')
title('Loop Input and Output - Full Range')
legend('Input', 'Output')

subplot (3,2,2)
plot(t(1:200)*1e3,Loop_input(1:200))
hold on
plot(t(1:200)*1e3,(Out_cos(1:200)), 'r')
ylabel('Amplitude')
xlabel('ms')
title('Loop Input and Output - Range 0 to 0.01 ms')
legend('Input', 'Output')

subplot (3,2,3)
plot(t(1601:1800)*1e3,Loop_input(1601:1800))
hold on
plot(t(1601:1800)*1e3,(Out_cos(1601:1800)), 'r')
ylabel('Amplitude')
xlabel('ms')
title('Loop Input and Output - Range 0.08 to 0.09 ms')
legend('Input', 'Output')

subplot (3,2,4)
plot(t(2001:2200)*1e3,Loop_input(2001:2200))
hold on
plot(t(2001:2200)*1e3,(Out_cos(2001:2200)), 'r')
ylabel('Amplitude')
xlabel('ms')
title('Loop Input and Output - Range 0.1 to 0.01 ms')
legend('Input', 'Output')

subplot (3,2,5)
plot(t(2701:2900)*1e3,Loop_input(2701:2900))
hold on
plot(t(2701:2900)*1e3,(Out_cos(2701:2900)), 'r')
ylabel('Amplitude')
xlabel('ms')
title('Loop Input and Output - Range 0.135 to 0.145 ms')
legend('Input', 'Output')

subplot (3,2,6)
plot(t(6001:6200)*1e3,Loop_input(6001:6200))
hold on
plot(t(6001:6200)*1e3,(Out_cos(6001:6200)), 'r')
ylabel('Amplitude')
xlabel('ms')
title('Loop Input and Output - Range 0.3 to 0.31 ms')
legend('Input', 'Output')

% Scope display

scope = dsp.TimeScope(...
        'Title', 'Loop Input(1) and Output(2)', ...
        'NumInputPorts', 2, ...
        'SampleRate', 20e+6, ...
        'ShowGrid', 1, ...
        'ShowLegend', 1, ...
        'BufferLength', 1e5, ...
        'TimeSpanOverrunAction', 'Wrap', ...
        'TimeSpan', 5e-6, ...
        'TimeUnits', 'Metric', ...
        'YLimits', [-2 2]);
    
for k=1:N           
           
    scope(Loop_input(k),Out_cos(k))
    for L=1:0.5e+6              % Delay for proper display
    end
    
end  
