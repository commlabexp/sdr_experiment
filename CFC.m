%% Coarse Frequency Correction

z = y.^Mod;       % raising the signal to its modulation
                  % order

pxx = pwelch(z,[],[],[],samplerate,'centered');  

figure(1)
pwelch(z,[],[],[],samplerate,'centered')    % Spectrum of the captured signal^Mod


maximum = max(pxx);                         % Maximum of pxx
szdim = length(pxx)                         % No of FFT bins
Bin_Resolution = samplerate/szdim           % Bin resolution (Hz)
Max_index=find(pxx==maximum)                % Index of the max FFT bin

CFoffset = ((Max_index - szdim/2) * (samplerate/szdim))/Mod       % CFoffset   

t=([0:length(y)-1]).'/samplerate;
corrected_y                                 % Captured y corrected by the coarse frequency error
        
r = y_CFC.^Mod;

figure(2) 
pwelch(r,[],[],[],samplerate,'centered')    % Spectrum of the corected signal^Mod


refConst = qammod(0:Mod-1,Mod,'UnitAveragePower',true);
cdCFC = comm.ConstellationDiagram('ReferenceConstellation',refConst, ...
    'SamplesPerSymbol',spsreceived,'SymbolsToDisplaySource','Property','SymbolsToDisplay',10000,'Title','Signal Constellation after CFC');

My_CFC = mean(abs(y_CFC));   % To normalize the constellatoon diagram

y_CFC = (y_CFC./My_CFC);

cdCFC(y_CFC)


