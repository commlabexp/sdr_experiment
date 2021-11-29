%% Symbol Synchronizer Loop
% ---------------------
% Implements Symbol Timing Recovery with configurable Timing Error Detector
% (TED) and Interpolator, using both a Proportional-plus-integrator (PI)
% Controller and a Modulo-1 Counter to control the interpolator. The TED
% can be configured as a Maximum-likelihood (ML) TED (MLTED), or Early-Late TED (ELTED), whereas the interpolator is a Parabolic
% Interpolator.

%% Parameters
L        = 4;          % Oversampling factor
M        = 4;          % Constellation order
nSymbols = 75000;
Bn_Ts    = 0.001;      % PLL noise bandwidth (Bn) times symbol period (Ts)
eta      = 1;          % PLL Damping Factor
rollOff  = 0.5;        % Pulse shaping roll-off factor
rcDelay  = 20;         % Receiver Square Raised cosine filter span
Ex       = 1;          % Average symbol energy
TED      = 'ELTED';    % TED Type (choose 'MLTED' , 'ELTED' or 'ZCTED')

%% Matched Filter (MF)
RXFILT  = comm.RaisedCosineReceiveFilter( ...
    'InputSamplesPerSymbol', L, ...
    'DecimationFactor',1, ...
    'RolloffFactor', rollOff, ...
    'FilterSpanInSymbols', rcDelay);    % Generating the square-root, raised cosine filter coefficients
        
mf  = RXFILT.coeffs.Numerator;

%% dMF

h = [0.5 0 -0.5]; % kernel function
central_diff_mf = conv(h, mf);
% Skip the filter delay
dmf = central_diff_mf(2:1+length(mf));

%% Signal at output of MF and dMF

% MF
rxSample = step(RXFILT,FFC_out);

%dMF
rxSampleDiff = filter(dmf, 1, FFC_out);


const  = qammod(0:M-1,M);           % Constellation
Ksym   = modnorm(const, 'avpow', Ex);
 
switch (TED)
            case 'MLTED' % Maximum Likelihood TED
                
                [ xx ] = MLTED_symTimingLoop(L, Ex, rxSample, rxSampleDiff, eta, Bn_Ts, ...
                       const);
                   
                             
            case 'ELTED' % Early Late TED
                [ xx ] = ELTED_symTimingLoop(L, Ex, rxSample, eta, Bn_Ts, ...
                       const);
        
  
end


cdy_Symbol = comm.ConstellationDiagram('ReferenceConstellation',refConst, ...
    'SamplesPerSymbol',1,'SymbolsToDisplaySource','Property','SymbolsToDisplay',10000, ...
    'Title','Signal after Symbol sync','EnableMeasurements',true);

Mxx = mean(abs(xx(end-round(0.1*nSymbols):end)));   % To normalize the constellatoon diagram

cdy_Symbol((xx(end-round(0.1*nSymbols):end)./Mxx))

%% Eye Diagram

% edy_filt = comm.EyeDiagram('SampleRate',samplerate*spsreceived,'SamplesPerSymbol',spsreceived,'SampleOffset',2*spsreceived);  % Eye diagram of the filtered signal
edy_filt = comm.EyeDiagram('SampleRate',samplerate,'SamplesPerSymbol',1);
edy_filt(xx./Mxx)
