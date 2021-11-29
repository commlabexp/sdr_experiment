% Capture Rx signal of PLUTO

%% Parameters

deviceName = 'Pluto';
samplerate = 40e6;              % Sample rate in PLUTO of the received signal
spsreceived = 4;                % Samples per received symbol
Frequency = 1.8e9;              % Frequency of the received signal
Mod = 4;                        % Modulation order
alpha = 0.5;                     % Rolloff factor
span = 20;                      % Transmit filter span in symbols
%% Capture
rx = sdrrx(deviceName,'BasebandSampleRate',samplerate, ...
    'CenterFrequency',Frequency,'OutputDataType','double');

capture(rx,3e5,'Filename','PlutoTxRx_alfa_05_SPS_4_10MBs_RRCspan_20_3x100K_bits_test.bb');     % The captured Rx signal is stored in rx_xxxx.bb

% release(rx);
 
% RadioSettings = info(rx)
