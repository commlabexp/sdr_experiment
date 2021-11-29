
% Test Signal

%% Parameters
Frequency = 1.8e9;       % Frequency of the received signal
spsreceived = 4;         % Samples per Symbol (sps)
Mod        = 4;          % Constellation order
nSymbols = 75000;        % Number of transmit symbols
samplerate = 40e6;       % Test signal sample rate (Hz)
PhaseOffset = 0;         % Phase offset of the test signal (Deg)
FrequencyOffset = 27500; % Frequency offset of the test signal (Hz)
%FrequencyOffset = 0;
%samplerateOffset = (5/3)*10^-5;   % 4.165 PPM (To match the Farrow algoritm for 300000 samples)
samplerateOffset = 0;
timingErr = 1;           % Delay (in samples) added
alpha  = 0.5;            % Pulse shaping roll-off factor
span  = 20;              % Transmit Square Raised cosine filter span
SNR   = 60;              % Test signal SNR (dB)

%%
% Tx Filter
TXFILT  = comm.RaisedCosineTransmitFilter( ...
    'OutputSamplesPerSymbol', spsreceived, ...
    'RolloffFactor', alpha, ...
    'FilterSpanInSymbols', span);

data = randi([0 Mod-1],nSymbols,1);       % Random QPSK symbols
modSig = pskmod(data,Mod,pi/4);           % QPSK modulation

txSig = TXFILT(modSig);
txSig = awgn(txSig,SNR-10*log10(spsreceived),'measured');         % Adding AWGN to the signal

fixedDelay = dsp.Delay(timingErr);

delaySig = fixedDelay(txSig);           % Delayed signal

newsamplerate = samplerate*(1+samplerateOffset);
fs1 = samplerate;
fs2 = newsamplerate;
LagrangeOrder = 2; % 2 = parabolic interpolation
frc = dsp.FarrowRateConverter('InputSampleRate',fs1,...
                              'OutputSampleRate',fs2,...
                              'PolynomialOrder',LagrangeOrder); % Sample rate converter

pfo = comm.PhaseFrequencyOffset('PhaseOffset',PhaseOffset, ...
    'FrequencyOffset',FrequencyOffset,'SampleRate',samplerate);

modSigOffset = pfo(delaySig);

newsamplerateSig = frc(modSigOffset);

y = newsamplerateSig;

%%
refConst = qammod(0:Mod-1,Mod,'UnitAveragePower',true);
cdTransmited = comm.ConstellationDiagram('ReferenceConstellation',refConst, ...
    'SamplesPerSymbol',spsreceived,'SymbolsToDisplaySource','Property','SymbolsToDisplay',10000,'Title','Test Signal Constellation');

My = mean(abs(y));   % To normalize the constellatoon diagram

y = (y./My);

cdTransmited(y)

