%%
%-------------------------Simulation Parameters--------------------
numberOfBits = log2(M);
blockLength = numberOfSymbols*numberOfBits; % o que é um bloco?

SNR = 30:-5:0;
modulationIndexVector = [0.001 0.01 0.05];

VDC = 3.25; 
maxAbsoluteValueModulation = 3;
maxModulationIndex = (maxLEDVoltage - VDC)/VDC;
%-------------------------Simulation Parameters--------------------
%%
%-------------------------LED Parameters-----------------------------------
fs = 2e6; 
Poptical = @(ledLuminousEfficacy,electricalPower,kNLin) (ledLuminousEfficacy.*electricalPower)./((1 + (ledLuminousEfficacy.*electricalPower./(maxLuminousIntensityLED/1000)).^(2*kNLin)).^(1/(2*kNLin)));
%-------------------------LED Parameters-----------------------------------
%-------------------------Photodiode Parameters----------------------------
A = 1e-4;   %photodiode area (cm)
R = 0.5;
%-------------------------Photodiode Parameters----------------------------

%-------------------------Pre Amplifier Parameters-------------------------
transimpedanceGain = 10; %not used
%-------------------------Pre Amplifier Parameters-------------------------
%%
%-------------------------Transmission Parameters--------------------------
if M == 4
    NumTx = 4; NumRx = 4; 
end
if M == 8
    NumTx = 8; NumRx = 8; 
end

kNonLinearity = 2;
H_0 = zeros(NumTx, NumRx);
LEDtxAngle = zeros (NumTx, NumRx); %deg2rad(0) = 0
PDrxAngle = zeros (NumTx, NumRx);  %deg2rad(0) = 0


% d =25e-2;  %distance between LED and photodiode (cm), altura de uma mesa
% R_rx = 10e-2; %fixed distance on receivers

%if ((LED_POSITION == 21)||(LED_POSITION == 22)||(LED_POSITION == 23)), d = 50e-2; end 
%if ((LED_POSITION == 21)||(LED_POSITION == 22)||(LED_POSITION == 23)), R_rx = R_tx; end 


if (LED_GEOMETRY == 1), R_rx = R_tx; end 


halfAngleLED = deg2rad(25);
FOV = deg2rad(25);
LambertianMode = (log(1/2))/log(cos(halfAngleLED));

if  (LED_GEOMETRY == 0)  % Circular Geometry
    alpha = zeros((NumTx/2) - 1);
    poligonSide = zeros((NumTx/2) - 1);
    for p = 1:(NumTx/2 - 1)
        alpha(p) = ((2*pi)/NumTx)*p; % p = 1,...,NumTx/2
        poligonSide(p) = sqrt(R_rx^2+R_tx^2+(2*R_rx*R_tx*cos(alpha(p))));
    end;
    for idxN = 2:NumRx
        if idxN <  ((NumRx/2)+1), LEDtxAngle(1, idxN) = atan(poligonSide(idxN-1)/d); end
        if idxN == ((NumRx/2)+1), LEDtxAngle(1, idxN) = atan((2*R_rx)/d); end
        if idxN >  ((NumRx/2)+1), LEDtxAngle(1, idxN) = atan(poligonSide(NumRx - (idxN-1))/d); end
    end
    LEDtxAngle(1, 1) = atan((R_tx-R_rx)/d); 
end

       
if  (LED_GEOMETRY == 0) && (NumTx == 4) && (NumRx == 4) % Circular Geometry
    LEDtxAngle = [LEDtxAngle(1, 1) LEDtxAngle(1, 2) LEDtxAngle(1, 3) LEDtxAngle(1, 4)
                  LEDtxAngle(1, 4) LEDtxAngle(1, 1) LEDtxAngle(1, 2) LEDtxAngle(1, 3);
                  LEDtxAngle(1, 3) LEDtxAngle(1, 4) LEDtxAngle(1, 1) LEDtxAngle(1, 2);
                  LEDtxAngle(1, 2) LEDtxAngle(1, 3) LEDtxAngle(1, 4) LEDtxAngle(1, 1)];
    PDrxAngle = LEDtxAngle;
end

if  (LED_GEOMETRY == 0) && (NumTx == 8) && (NumRx == 8) % Circular Geometry
    LEDtxAngle = [LEDtxAngle(1, 1) LEDtxAngle(1, 2) LEDtxAngle(1, 3) LEDtxAngle(1, 4) LEDtxAngle(1, 5) LEDtxAngle(1, 6) LEDtxAngle(1, 7) LEDtxAngle(1, 8)
                  LEDtxAngle(1, 8) LEDtxAngle(1, 1) LEDtxAngle(1, 2) LEDtxAngle(1, 3) LEDtxAngle(1, 4) LEDtxAngle(1, 5) LEDtxAngle(1, 6) LEDtxAngle(1, 7);
                  LEDtxAngle(1, 7) LEDtxAngle(1, 8) LEDtxAngle(1, 1) LEDtxAngle(1, 2) LEDtxAngle(1, 3) LEDtxAngle(1, 4) LEDtxAngle(1, 5) LEDtxAngle(1, 6);                
                  LEDtxAngle(1, 6) LEDtxAngle(1, 7) LEDtxAngle(1, 8) LEDtxAngle(1, 1) LEDtxAngle(1, 2) LEDtxAngle(1, 3) LEDtxAngle(1, 4) LEDtxAngle(1, 5);                 
                  LEDtxAngle(1, 5) LEDtxAngle(1, 6) LEDtxAngle(1, 7) LEDtxAngle(1, 8) LEDtxAngle(1, 1) LEDtxAngle(1, 2) LEDtxAngle(1, 3) LEDtxAngle(1, 4);                 
                  LEDtxAngle(1, 4) LEDtxAngle(1, 5) LEDtxAngle(1, 6) LEDtxAngle(1, 7) LEDtxAngle(1, 8) LEDtxAngle(1, 1) LEDtxAngle(1, 2) LEDtxAngle(1, 3);                 
                  LEDtxAngle(1, 3) LEDtxAngle(1, 4) LEDtxAngle(1, 5) LEDtxAngle(1, 6) LEDtxAngle(1, 7) LEDtxAngle(1, 8) LEDtxAngle(1, 1) LEDtxAngle(1, 2);                  
                  LEDtxAngle(1, 2) LEDtxAngle(1, 3) LEDtxAngle(1, 4) LEDtxAngle(1, 5) LEDtxAngle(1, 6) LEDtxAngle(1, 7) LEDtxAngle(1, 8) LEDtxAngle(1, 1)];                  
      PDrxAngle = LEDtxAngle;
end
 

if (LED_GEOMETRY == 1) && (NumTx == NumRx) % Em linha
    for idxM = 1:NumTx
        for idxN = 1:NumRx
            if idxM ~= idxN
                LEDtxAngle(idxM, idxN) = atan((R_tx*abs(idxN-idxM))/d);   
                PDrxAngle(idxM, idxN)  = atan((R_tx*abs(idxM-idxN))/d);
            end
            if idxM == idxN
                LEDtxAngle(idxM, idxN) = 0;
                PDrxAngle(idxM, idxN)  = 0;
            end
        end
    end
end

for m = 1:NumTx
    for n = 1:NumRx
        H_0(m, n) = (LambertianMode + 1)/(2*pi) * (A/d^2) * (cos(LEDtxAngle(m,n)))^LambertianMode * cos(PDrxAngle(m, n)) * rectangularPulse(-1,1,LEDtxAngle(m, n)/FOV);             
    end
end

%-------------------------Transmission Parameters--------------------------
%%
%-------------------------Defining Variables Size--------------------------
INCREMENT = 100; % utilizado para defifinir o tamanho de bitInputData

berAux = zeros(monteCarloLoops,1);
ber = zeros(length(SNR),length(modulationIndexVector));

berSymbol = zeros(length(SNR),length(modulationIndexVector));
berPosition = zeros(length(SNR),length(modulationIndexVector));
berAuxSymbol = zeros(monteCarloLoops,1);
berAuxPosition = zeros(monteCarloLoops,1);


berInfo = zeros(length(SNR),length(modulationIndexVector));
berAuxTraining = zeros(monteCarloLoops,1);
berTraining = zeros(length(SNR),length(modulationIndexVector));


%pilotVector = zeros((blockLength+INCREMENT)/4, NumTx);
blockLengthFilter = ((blockLength+INCREMENT)*numberOfBits);

pulseShapingFilter = zeros(blockLengthFilter, NumTx);
pilotMatchedFilter = zeros((blockLength+INCREMENT)/2, NumTx);

transImpedanceGain = zeros(1, NumTx);
decDemodSignal = zeros((blockLength+INCREMENT)/2, NumRx);
bitOutputDataMO = zeros((blockLength+INCREMENT)/2, 2, NumTx);
%--------------------------------------------------------------------------
bitInputDataLength = (ceil((blockLength + INCREMENT)/8) * 4)/2;
convLength = bitInputDataLength + 1000 -1;         % 1000??
NFFT = 2^nextpow2(convLength);
%--------------------------------------------------------------------------
LEDRespVector = zeros(NFFT, NumTx);
filteredVinAux = zeros(NFFT, NumTx);
filteredVin = zeros(bitInputDataLength, NumTx);
iLEDOutput = zeros(bitInputDataLength, NumTx);
eletricalPowerOutput = zeros(bitInputDataLength, NumTx);
opticalPowerOutput = zeros(bitInputDataLength, NumTx);
opticalPowerInput = zeros(bitInputDataLength, NumRx);
receivedVoltageSignal = zeros (bitInputDataLength, NumRx);
receivedCurrentSignalAC = zeros(bitInputDataLength, NumRx);
receivedCurrentSignalPower = zeros(NumRx, NumTx);
equalizedSignal = zeros(bitInputDataLength, 1);


%pulseShapingFilter = zeros(length(Vin), NumTx);
nAux = zeros(bitInputDataLength, NumRx);
nVector = zeros(bitInputDataLength, NumRx);
%yArray = zeros(blockLength, NumTx);

encodedSymbol = zeros(bitInputDataLength , NumTx);
decodedPosition = zeros(bitInputDataLength , numberOfBits);

%-------------------------Defining Variables Size--------------------------
%%
%-------------------------Debug--------------------------
WHITE_NOISE             = 1;% 
CHANNEL_EFFECT          = 1;% 
LED_BYPASS              = 0;%
%-------------------------Debug--------------------------
%%

