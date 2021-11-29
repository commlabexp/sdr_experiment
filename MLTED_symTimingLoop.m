function [ xx ] = MLTED_symTimingLoop(L, Ex, mfOut, dMfOut, eta, Bn_Ts, const)

%
% Input Arguments:
% L       -> Oversampling Factor (SpS)
% Ex      -> Average symbol energy
% mfOut   -> MF output sequence sampled at L samples/symbol
% dMfOut  -> Derivative MF output sequence sampled at L samples/symbol
%            (used only for the Maximum-Likelihood TED)
% eta     -> PLL Damping Factor
% Bn_Ts   -> PLL noise bandwidth (Bn) times symbol period (Ts)
% const   -> Modulation Constellation


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

Kp = 0.23;
K  = 1;         % Assume channel gain is unitary
Kp = K*Ex*Kp;

% Counter Gain
K0 = -1;
% Note: this is analogous to VCO or DDS gain, but in the context of timing
% recovery loop.

% PI Controller Gains:
[ K1, K2 ] = timingLoopPIConstants(Kp, K0, eta, Bn_Ts, L)

%% Timing Recovery Loop

% Constants
nSamples = length(mfOut);
nSymbols = ceil(nSamples / L);

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
W1        = 1/L;
FX1 =0; FX2 =0; FX3 =0; FX4 =0; FX5 =0;
FDX1 =0; FDX2 =0; FDX3 =0; FDX4 =0; FDX5 =0;

for n=1:length(mfOut)
 
        
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
                    tempFx = -0.5*mfOut(n);
                    tempFdx = -0.5*dMfOut(n);
                    
                    % Interpolants (See Eq. 8.77)
                    vF2 = -tempFx + FX1 +FX2 - FX3;
                    vF1 = tempFx - FX1 + FX4 + FX2 + FX3;
                    vF0 = FX5;
                    xi = (vF2*mu + vF1)*mu + vF0;
                    
                    vF2 = -tempFdx + FDX1 +FDX2 - FDX3;
                    vF1 = tempFdx - FDX1 + FDX4 + FDX2 + FDX3;
                    vF0 = FDX5;
                    dxi = (vF2*mu + vF1)*mu + vF0;
    
                    %% Timing Error Detector:
                %e = real(xi) * real(dxi) + ...
                    %imag(xi) * imag(dxi);
                e = sign(real(xi)) * real(dxi) + ...
                    sign(imag(xi)) * imag(dxi);
                
                % Save interpolant for MF output, Timing Error and Fractional
        % Interval for plotting:
        xx(k)   = xi;
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
    w    = 1/L + v;      % Adjust Counter Step
    
    %%% Compute next states
    FX3 = FX2;
    FX2 = FX1;
    FX1 = -0.5*mfOut(n);
    FX5 = FX4;
    FX4 = mfOut(n);
    
    FDX3 = FDX2;
    FDX2 = FDX1;
    FDX1 = -0.5*dMfOut(n);
    FDX5 = FDX4;
    FDX4 = dMfOut(n);
    
    V21 = v2;
    
    NCO1 = nco;
    MU1 = mu;
    W1 = w;
            
   
end
%% Dynamic Debug Plots

% Debug using scope
for k=3:nSymbols*0.25           % Display only the leading 25% of the trace (to save display time)
    
                hTScopeCounter(mu_k(k));
       
                hScope(xx(k)/1.47)
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
