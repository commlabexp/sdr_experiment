%% Fine Carrier Synchronizer

% General system details 
DampingFactor = 0.9;
Bn = 0.001;                 % Normalized to the symbol rate (10MHz*0.001=10KHz)
K0 = -1;                    % DDS Gain (Negative feedback)
Kp = 2;                     % Phase Detector Gain
delta_time = 1/samplerate;

 % Real and Imaginary parts of the input signal y_filt
 I = real(y_CFC);
 Q = imag(y_CFC);
 
 % DDS definitions
 theta = 0;             % Initial sinusoid arguments
 DDSout = 0;            % DDS integrator initial value
 vi = 0;                % Integrator initial value
 
 % Calculate coefficients for FFC loop filter

% K1
K1 = (4*DampingFactor*Bn)/(Kp*(spsreceived*(DampingFactor+1/(4*DampingFactor))));

% K2
K2 = (4*Bn^2)/(Kp*(spsreceived^2)*((DampingFactor+1/(4*DampingFactor))^2));

 % Allocate room for vectors used by the subsystem
 Iprim = zeros(1, length(I));
 Qprim = zeros(1, length(Q));
 d0 = zeros(1, length(I));
 d1 = zeros(1, length(Q));
 e = zeros(1, length(I));
 thetaplot = zeros(1, length(I));
 

 for k = 1:length(I)
 % Counter Clockwise Phase Rotation
 Iprim(k) = cos(theta)*I(k)+sin(theta)*Q(k);
 Qprim(k) = -sin(theta)*I(k)+cos(theta)*Q(k);

 % QPSK symbol decisions
 d0(k) = sign(Iprim(k));
 d1(k) = sign(Qprim(k));

 % Phase error detector
 e(k) = Qprim(k)*d0(k)-Iprim(k)*d1(k);

 % Proportional Plus Integrator loop filter
 filtIn = Kp*e(k);      % Filter input
 vp = K1*filtIn;        % Proportional component
 vi = vi+K2*filtIn;     % Integrator component
 filtOut = vp+vi;       % Filter output

% Direct Digital Synthesizer
 DDSin = K0*filtOut;        % DDS input
 DDSout = DDSout+DDSin;     % DDS ouput = old DDS value + new
 theta = DDSout;            % New phase compensation

 thetaplot(k) = theta;
  end
 
 figure (3)
 plot(thetaplot)
 title('DDS Phase Output')
 
 Frequency_correction = (thetaplot(end)-thetaplot(end-100000))/(delta_time*100000*2*pi)
 
 FFC_out = transpose(complex(Iprim,Qprim));     % Transpose column vector to raw vector
 FFC_out = FFC_out*exp(1i*pi/4);                % Shift by 45 deg
  
cdy_FFC = comm.ConstellationDiagram('ReferenceConstellation',refConst, ...
    'SamplesPerSymbol',spsreceived,'SymbolsToDisplaySource','Property','SymbolsToDisplay',10000,'Title','Signal after FFC');

MFFC_out = mean(abs(FFC_out));   % To normalize the constellatoon diagram

FFC_out = (FFC_out./MFFC_out);

cdy_FFC(FFC_out)


%% Eye Diagram

%edy_filt = comm.EyeDiagram('SampleRate',samplerate,'SamplesPerSymbol',spsreceived);
%edy_filt(FFC_out)