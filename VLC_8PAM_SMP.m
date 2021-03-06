%   MIMO implementation with Spaial multiplexing (SMP) technique, using PAM symbols
%   Spatial multiplexing - different signal from each transmitter

%%
%-------------------------Defining Variables Size--------------------------
INCREMENT = 100; % utilizado para defifinir o tamanho de binaryInputData

berAux = zeros(monteCarloLoops, NumRx);
ber = zeros(length(SNR),length(modulationIndexVector));

berAuxInfo = zeros(monteCarloLoops,1);
berInfo = zeros(length(SNR),length(modulationIndexVector));
berAuxTraining = zeros(monteCarloLoops,1);
berTraining = zeros(length(SNR),length(modulationIndexVector));


pilotVector = zeros((blockLength+INCREMENT)/2, NumTx);

transImpedanceGain = zeros(1, NumTx);
equalizedSignal = zeros((blockLength+INCREMENT)/2, NumTx);
decDemodSignal = zeros((blockLength+INCREMENT)/2, NumTx);
binaryInputDataMO = zeros((blockLength+INCREMENT)/2, 3, NumTx);
binaryOutputDataMO = zeros((blockLength+INCREMENT)/2, 3, NumTx);
%--------------------------------------------------------------------------
binaryInputData = randi([0,1],(blockLength + INCREMENT)/2,3); 
%binaryInputData = reshape(binaryInputData,[],numberOfBits);
deciInputData = bi2de(binaryInputData); 
Vin = deciInputData; 
convLength = length(Vin) + 1000 -1;         % 1000??
NFFT = 2^nextpow2(convLength);
%--------------------------------------------------------------------------
LEDRespVector = zeros(NFFT, NumTx);
filteredVinAux = zeros(NFFT, NumTx);
filteredVin = zeros(length(Vin), NumTx);
iLEDOutput = zeros(length(Vin), NumTx);
eletricalPowerOutput = zeros(length(Vin), NumTx);
opticalPowerOutput = zeros(length(Vin), NumTx);
opticalPowerInput = zeros(length(Vin), NumRx);
receivedVoltageSignal = zeros (length(Vin), NumRx);
receivedCurrentSignalAC = zeros(length(Vin), NumRx);
receivedCurrentSignalPower = zeros(NumTx, NumRx);

%pulseShapingFilter = zeros(length(Vin), NumTx);
nAux = zeros(length(Vin), NumRx);
nVector = zeros(length(Vin), NumRx);
%yArray = zeros(blockLength, NumTx);
%-------------------------Defining Variables Size--------------------------
 

%%
for index = 1:length(modulationIndexVector) 
    modulationIndex = modulationIndexVector(index);
    if modulationIndex > maxModulationIndex
        warning('Modulation Index may cause undesired nonlinear effects')
    end
    maxVoltage = VDC*(1+modulationIndex);
    deltaV = maxVoltage - VDC;
    VoltageConstant = modulationIndex*maxVoltage/((1+modulationIndex)*maxAbsoluteValueModulation); 
    for SNRIndex = 1:length(SNR)
        for j = 1:monteCarloLoops
            j;
            for k = 1:NumTx
                binaryInputData = randi([0,1],(blockLength + INCREMENT)/2,3); 
                %binaryInputData = reshape(binaryInputData,[],numberOfBits);
                
                binaryInputDataMO(:,:,k) = binaryInputData;
                deciInputData = bi2de(binaryInputData); 

                pilot = pammod(deciInputData,2^numberOfBits,0,'gray');
            
                pilotVector(:,k) = pilot; 
            end
            Vin = pilotVector;

            convLength = length(Vin) + 1000 -1;         % 1000??
            NFFT = 2^nextpow2(convLength);
            VinFreq = fft(Vin,NFFT);                    %Colocando o Sinal no dom�nio da frequ�ncia         
            fR = fs/2 * linspace(0,1,NFFT/2 + 1)*2*pi;
            wR = [-fliplr(fR(2:end-1)) fR];
            LEDResp = freqRespLED(wR);                  %Resposta em frequ�ncia do LED

            for k = 1:NumTx
                LEDRespVector(:,k) = LEDResp;
                filteredVinAux(:,k) = real(ifft(VinFreq(:,k).*fftshift(LEDRespVector(:,k))));
                filteredVin(:,k) = filteredVinAux(1:length(Vin),k); 
                VoltageConstant(:,k) = modulationIndex.*maxVoltage/((1+modulationIndex)*max(filteredVin(:,k)));             
                filteredVin(:,k) = (filteredVin(:,k).*VoltageConstant(:,k)) + VDC;           
                iLEDOutput(:,k) = I_V_Fun(filteredVin(:,k),VT,nLED,ISat);              
                eletricalPowerOutput(:,k) = filteredVin(:,k).*iLEDOutput(:,k);   
                opticalPowerOutput(:,k) = Poptical(ledLuminousEfficacy,eletricalPowerOutput(:,k),kNonLinearity);
            end
            
            %opticalPowerOutputConvolved = opticalPowerOutput;
            
%-------------------------LED Bypass--------------------------
%%
            if LED_BYPASS == 1
                opticalPowerOutput = Vin;
            end
%-------------------------LED Bypass--------------------------
%%
%-------------------------Channel Effect--------------------------
%%
            if CHANNEL_EFFECT == 1
                opticalPowerInput = opticalPowerOutput*H_0;
            end
%-------------------------Channel Effect--------------------------
%%
            receivedCurrentSignal = opticalPowerInput*R*A;

            for k = 1:NumRx
                receivedCurrentSignalAC(:,k) = receivedCurrentSignal(:,k) - mean(receivedCurrentSignal(:,k));
            end
            
            receivedCurrentSignalPower = sum(receivedCurrentSignalAC.*receivedCurrentSignalAC, 1)/length(receivedCurrentSignal); % como isso da pot�ncia?

            %receivedCurrentSignalPower = receivedCurrentSignalPower';
            receivedVoltageSignalAux = receivedCurrentSignal; 

%-------------------------White Noise--------------------------
%%
            if WHITE_NOISE == 1
                n = randn(length(opticalPowerOutput),1); 
                powerNoiseAux = n'*n/(length(n));
                powerNoise = (receivedCurrentSignalPower/db2pow(SNR(SNRIndex)));
                PowerNoise = powerNoise/powerNoiseAux;
                for k = 1:NumRx
                    PowerNoiseAux = PowerNoise(1,k);
                    nAux(:,k) = n.*sqrt(PowerNoiseAux);               
                    nVector(:,k) = nAux(:,k);
                end
                receivedVoltageSignalAux = (receivedCurrentSignal + nVector);
            end
%-------------------------White Noise--------------------------
%%                    
            pilotVariance = sum(var(pilotVector))/size(pilotVector,2);
            for k = 1:NumRx                
                receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux(:,k));            
                transImpedanceGain(:,k) = maxAbsoluteValueModulation/max(receivedVoltageSignalAux(:,k));            
                %transImpedanceGain pra que serve?
                debug = sqrt(pilotVariance/var(receivedVoltageSignalAux(:,k)));
                %receivedVoltageSignal(:,k) =  receivedVoltageSignalAux(:,k)*sqrt(pilotVariance/var(receivedVoltageSignalAux(:,k)));
                receivedVoltageSignal(:,k) =  receivedVoltageSignalAux(:,k)*debug;
            end
            output = receivedVoltageSignal;
%-------------------------Least-Squares Equalization--------------------------
%%          
            unbiasedReceivedVoltageSignal = receivedVoltageSignal - VDC;                   %pra que serve esse sinal?? 
            delay = 0; 
            
            if LS_EQUALIZATION == 1 % equaliza��o COM sequ�ncia de treinamento
            	[f] = LSEqualizer(index, pilot, receivedVoltageSignal, numberOfBits, COD_GEOMETRY, LED_POSITION);                              
                output = filter(f, 1, receivedVoltageSignal) ;  
                                     

                %auxEqualizerLength(j) = equalizerLength; 
            end      
%-------------------------Least-Squares Equalization--------------------------
%%

%-------------------------Demodulation--------------------------
%%     
            for k = 1:NumRx                
                equalizedSignal(:,k) = pam8HardThreshold(output(:,k));
            
                [corr,lags] = xcorr(equalizedSignal(:,k), receivedVoltageSignal(:,k));
                [~,idx] = max(abs(corr));
                delay(:,k) = abs(lags(idx));
            
                decDemodSignal(:,k) = pamdemod(equalizedSignal(:,k),2^numberOfBits,0,'gray'); %sinal demodulado               
                binaryOutputDataMO(:,:,k) = de2bi(decDemodSignal(:,k) ,numberOfBits); 
            
%-------------------------Demodulation--------------------------

                berAux(j, k) = sum(sum(abs(binaryOutputDataMO(delay(:,k)+1:(blockLength/2) + delay(:,k),:, k) - binaryInputDataMO(1:(blockLength/2),:,k))))./(blockLength*NumRx);
            end
            
        end
        ber(SNRIndex,index) = mean(mean(berAux));

    end
end

save(['.' filesep 'resultsBER' filesep FileName '.mat'],'ber','SNR', 'monteCarloLoops', 'numberOfSymbols', 'modulationIndexVector');