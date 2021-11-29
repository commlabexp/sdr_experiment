% Read thee captured Rx signal of PLUTO

%% Parameters

deviceName = 'Pluto';
samplerate = 40e6;              % Sample rate in PLUTO of the received signal
spsreceived = 4;                % Samples per received symbol
Frequency = 1.8e9;              % Frequency of the received signal
Mod = 4;                        % Modulation order
alpha = 0.5;                    % Rolloff factor
span = 20;                      % Transmit filter span in symbols

%% Read the captured data

bbr = comm.BasebandFileReader('PlutoTxRx_alfa_05_SPS_4_10MBs_RRCspan_20_3x100K_bits_test.bb');         % Read the captured file
% info(bbr)                                                 % More info

y = [];

while ~isDone(bbr)                                          % Concatenate all the samples of the captured data as vector y
    x = bbr();
    y = cat(1,y,x);
end
% bbrinfo = bbr.info;                                       % Info if all
                                                            % the samples were read
release(bbr)                                                % Release the baseband file reader resources
%%

refConst = qammod(0:Mod-1,Mod,'UnitAveragePower',true);
cdReceived = comm.ConstellationDiagram('ReferenceConstellation',refConst, ...
    'SamplesPerSymbol',spsreceived,'SymbolsToDisplaySource','Property','SymbolsToDisplay',100000,'Title','Received Signal Constellation');

y_max = max(abs(y));

cdReceived(y./y_max)

