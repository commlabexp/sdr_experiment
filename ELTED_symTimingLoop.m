function [ xx ] = ELTED_symTimingLoop(L, Ex, mfOut, eta, Bn_Ts, const)

%
% Input Arguments:
% L       -> Oversampling Factor
% Ex      -> Average symbol energy
% mfOut   -> MF output sequence sampled at L samples/symbol
% eta     -> PLL Damping Factor
% Bn_Ts   -> PLL noise bandwidth (Bn) times symbol period (Ts)
% const   -> Modulation Constellation

x = resample(mfOut,1,L/2);      % ELTED works best with 2 samples per symbol (sps).
Ld = 2;                         % Down sampled samples per symbol
%% System Objects for Step-by-step Debugging of the Loop

% Constellation Diagram

    hScope = comm.ConstellationDiagram(...
        'SymbolsToDisplaySource', 'Property',...
        'SamplesPerSymbol', 1, ...
        'Position' , [1150 325 580 425], ...
        'ShowReferenceConstellation', false);
        
    hScope.XLimits = [-2 2]*max(real(const));
    hScope.YLimits = [-2 2]*max(imag(const));

% Time scope used to debug the fractional error

    hTScopeCounter = dsp.TimeScope(...
        'Title', 'Fractional Inverval', ...
        'NumInputPorts', 1, ...
        'ShowGrid', 1, ...
        'ShowLegend', 1, ...
        'BufferLength', 1e5, ...
        'TimeSpanOverrunAction', 'Wrap', ...
        'TimeSpan', 2e4, ...
        'TimeUnits', 'None', ...
        'YLimits', [-1 1]);

%% PLL Design

Kp = 2.7;
K  = 1;         % Assume channel gain is unitary
Kp = K*Ex*Kp;

% Counter Gain
K0 = -1;
% Note: this is analogous to VCO or DDS gain, but in the context of timing
% recovery loop.

% PI Controller Gains:
[ K1, K2 ] = timingLoopPIConstants(Kp, K0, eta, Bn_Ts, Ld)

%% Timing Recovery Loop

% Constants
nSamples = length(x);
nSymbols = ceil(nSamples / Ld);

% Preallocate
xx   = zeros(nSymbols, 1);
ee   = zeros(nSymbols, 1);
mu_k = zeros(nSymbols, 1);
vv   = zeros(nSymbols, 1);

% Initialize
k         = 1;
MU1       = 0;
NCO1      = 1;
V21       = 0;
W1        = 1/Ld;
F0X1 =0; F0X2 =0; F0X3 =0; F0X4 =0; F0X5 =0;
F1X1 =0; F1X2 =0; F1X3 =0; F1X4 =0; F1X5 =0;
F2X1 =0; F2X2 =0; F2X3 =0; F2X4 =0; F2X5 =0;
X1 =0; X2 =0;

for n=3:length(x)
 
        
    %% Parabolic interpolator
                % Evaluate arithmetic expresions in topological order
                
                % Compute nco, strobe, mu
                temp = NCO1 - W1;
                if temp < 0
                    strobe =1;
                    mu = NCO1/W1;
                    nco = 1 + temp;
                else
                    strobe = 0;
                    mu = MU1;
                    nco = temp;
                end
                % if strobe, make decisions and compute timing error
                if strobe == 0
                    e = 0;
                else
                    % Compute interpolant xi from MF outputs
                    
                    tempFx = -0.5*x(n);
                    vF2 = -tempFx + F0X1 +F0X2 - F0X3;
                    vF1 = tempFx - F0X1 + F0X4 + F0X2 + F0X3;
                    vF0 = F0X5;
                    xi = (vF2*mu + vF1)*mu + vF0;
                    
                    % Compute interpolant xi1 from MF outputs
                    tempFx = -0.5*X1;
                    vF2 = -tempFx + F1X1 +F1X2 - F1X3;
                    vF1 = tempFx - F1X1 + F1X4 + F1X2 + F1X3;
                    vF0 = F1X5;
                    xi1 = (vF2*mu + vF1)*mu + vF0;
    
                    % Compute interpolant xi2 from MF outputs
                    tempFx = -0.5*X2;
                    vF2 = -tempFx + F2X1 +F2X2 - F2X3;
                    vF1 = tempFx - F2X1 + F2X4 + F2X2 + F2X3;
                    vF0 = F2X5;
                    xi2 = (vF2*mu + vF1)*mu + vF0;
                    
                    %% Timing Error Detector:
                
                e = sign(real(xi1))*((real(xi))-(real(xi2))) + ...
                    sign(imag(xi1))*((imag(xi))-(imag(xi2)));
                
                % Save interpolant for MF output, Timing Error and Fractional
        % Interval for plotting:
        xx(k)   = xi1;
        ee(k)   = e;
        mu_k(k) = mu;
        vv(k)   = v;

        % Update Interpolant Index
        k = k+1;
                
                end
        
        
        %% Loop Filter
    v1   = K1*e;       % Proportional
    v2   = V21 + K2*e;  % Integral
    v    = v1 + v2;    % PI Output

    %% Modulo-1 Counter
    
    % Compute NCO input
    w    = 1/Ld + v;      % Adjust Counter Step
    
    %%% Compute next states
    F0X3 = F0X2;
    F0X2 = F0X1;
    F0X1 = -0.5*x(n);
    F0X5 = F0X4;
    F0X4 = x(n);
    
    F1X3 = F1X2;
    F1X2 = F1X1;
    F1X1 = -0.5*X1;
    F1X5 = F1X4;
    F1X4 = X1;
    
    F2X3 = F2X2;
    F2X2 = F2X1;
    F2X1 = -0.5*X2;
    F2X5 = F2X4;
    F2X4 = X2;
    
    X2 = X1;
    X1 = x(n);
    
    V21 = v2;
    
    NCO1 = nco;
    MU1 = mu;
    W1 = w;
            
   
end
%% Dynamic Debug Plots

% Debug using scope
for k=3:nSymbols*0.25           % Display only the leading 25% of the trace (to save display time)
    
    
            hTScopeCounter(mu_k(k));
        
            hScope(xx(k)/1.54)
    
end   
%% Static Debug Plots

    figure
    plot(ee)
    ylabel('Timing Error $e(t)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')

    figure
    plot(vv)
    title('PI Controller Output')
    ylabel('$v(n)$', 'Interpreter', 'latex')
    xlabel('Sample $n$', 'Interpreter', 'latex')

    figure
    plot(mu_k, '.')
    title('Fractional Error')
    ylabel('$\mu(k)$', 'Interpreter', 'latex')
    xlabel('Symbol $k$', 'Interpreter', 'latex')


end
