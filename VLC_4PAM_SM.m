%   MIMO implementation with Spatial Modulation (SM) technique, using PAM symbols
%   Spatial Modulation - MIMO and digital modulation technique
 
%%
%-------------------------Pulse Shapping--------------------------
% (span * numSamplesPerSymbol) must be even
rolloff = 0.25;             % Rolloff factor
span = 6;                   % Filter span in symbols
numSamplesPerSymbol = 4;    % Samples per symbol

rCosineCoeff = rcosdesign(rolloff, span, numSamplesPerSymbol); %Generate the square-root, raised cosine filter coefficients.
%-------------------------Pulse Shapping--------------------------
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
blockLengthFilter = ((blockLength+INCREMENT)*numberOfBits) + 1 + (numSamplesPerSymbol/2)*10;

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

            lengthbitInputData = ceil((blockLength + INCREMENT)/8) * 8;
            bitInputData = randi([0,1],lengthbitInputData,1); % porque +100??

            %SM encoder, signal index & transmitter index

            bitInputData = reshape(bitInputData,[],numberOfBits);
            
            %symbolEncoder = zeros (length(bitInputData)/2, 2);
            %positionEncoder = zeros (length(bitInputData)/2, 2);
            
            symbolDecoder = zeros (length(bitInputData)/2, 2);
            positionDecoder = zeros (length(bitInputData)/2, 2);
            
            bitInputDataSymbol = zeros (length(bitInputData), 2);
            bitInputDataPosition = zeros (length(bitInputData), 2);
            
            bitInputDataSymbol(1:2:end, 1:2) = bitInputData(1:2:end, 1:2);
            bitInputDataPosition(2:2:end, 1:2) = bitInputData(2:2:end, 1:2);
            
            symbolEncoder = bitInputData(1:2:end, 1:2);
            positionEncoder = bitInputData(2:2:end, 1:2);

            deciInputData = bi2de(symbolEncoder); 
            pilot = pammod(deciInputData,2^numberOfBits,0,'gray');
           

            %Vin = pilot;
            %buffer = zeros (4, lengthVin);
            %Vin = [Vin; zeros(max(0,4*lengthVin -length(Vin)),1)];

%-------------------------Encoder--------------------------
%% 
            for k = 1:size(positionEncoder,1)
                switch positionEncoder (k,1)
                    case 0
                        switch positionEncoder (k, 2)
                            case 0
                                encodedSymbol(k, :) = [pilot(k);0;0;0];
                            case 1
                                encodedSymbol(k, :) = [0;pilot(k);0;0];
                        end
                    case 1
                        switch positionEncoder (k, 2)
                            case 1
                                encodedSymbol(k, :) = [0;0;pilot(k);0];
                            case 0
                                encodedSymbol(k, :) = [0;0;0;pilot(k)];
                        end
                end
            end
%-------------------------Encoder--------------------------
%%
            Vin = encodedSymbol;

            convLength = length(Vin) + 1000 -1;         % 1000??
            NFFT = 2^nextpow2(convLength);
            VinFreq = fft(Vin,NFFT);                    %Colocando o Sinal no domínio da frequência         
            fR = fs/2 * linspace(0,1,NFFT/2 + 1)*2*pi;
            wR = [-fliplr(fR(2:end-1)) fR];
            LEDResp = freqRespLED(wR);                  %Resposta em frequência do LED
            
            for count = 1:NumTx
                LEDRespVector(:,count) = LEDResp;
                filteredVinAux(:,count) = real(ifft(VinFreq(:,count).*fftshift(LEDRespVector(:,count))));
                filteredVin(:,count) = filteredVinAux(1:length(Vin),count); 
                VoltageConstant(:,count) = modulationIndex.*maxVoltage/((1+modulationIndex)*max(filteredVin(:,count)));             
                filteredVin(:,count) = (filteredVin(:,count).*VoltageConstant(:,count)) + VDC;  
                iLEDOutput(:,count) = I_V_Fun(filteredVin(:,count),VT,nLED,ISat);              
                eletricalPowerOutput(:,count) = filteredVin(:,count).*iLEDOutput(:,count);   
                opticalPowerOutput(:,count) = Poptical(ledLuminousEfficacy,eletricalPowerOutput(:,count),kNonLinearity);
            end
            
            %fazer um chaveamente que desliga os LEDs que não transmitem
            %informação

            switchOff = zeros(bitInputDataLength , NumTx);
            [row,col] = find (encodedSymbol);            
            for rowCount=1:(size(encodedSymbol,1))
                switchOff(row(rowCount),col(rowCount)) = 1;
            end
            
            %opticalPowerOutput = opticalPowerOutput .* switchOff;
            
            %opticalPowerOutputConvolved = opticalPowerOutput;
%-------------------------LED Bypass--------------------------
%%
            if LED_BYPASS == 1
                opticalPowerOutput = Vin;
            end
%-------------------------LED Bypass--------------------------
%%
            opticalPowerInput = opticalPowerOutput;
%-------------------------Channel Effect--------------------------
%%           
            if CHANNEL_EFFECT == 1
                opticalPowerInput = opticalPowerOutput*H_0;
            end
%-------------------------Channel Effect--------------------------
%%
            receivedCurrentSignal = opticalPowerInput*R*A;
            
%             for count = 1:NumRx
%                 receivedCurrentSignalAC(:,count) = receivedCurrentSignal(:,count) - mean(receivedCurrentSignal(:,count));
%                 receivedCurrentSignalPower(:,count) = (receivedCurrentSignalAC(:,count)'*receivedCurrentSignalAC(:,count) )/length(receivedCurrentSignal); % como isso da potência?
%             end
            
            
            for k = 1:NumRx
                receivedCurrentSignalAC(:,k) = receivedCurrentSignal(:,k) - mean(receivedCurrentSignal(:,k));
            end
            
            receivedCurrentSignalPower = sum(receivedCurrentSignalAC.*receivedCurrentSignalAC, 1)/length(receivedCurrentSignal); % como isso da potência?

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
            for count = 1:NumRx                
                receivedVoltageSignalAux = receivedVoltageSignalAux - mean(receivedVoltageSignalAux(:,count));            
                transImpedanceGain(:,count) = maxAbsoluteValueModulation/max(receivedVoltageSignalAux(:,count));            
                %transImpedanceGain pra que serve?
                debug = sqrt(var(pilot)/var(receivedVoltageSignalAux(:,count)));
                %receivedVoltageSignal(:,count) =  receivedVoltageSignalAux(:,count)*sqrt(var(pilot)/var(receivedVoltageSignalAux(:,count)));
                receivedVoltageSignal(:,count) =  receivedVoltageSignalAux(:,count)*debug;
            end
            %equalizedSignal = receivedVoltageSignal;          
%%

%-------------------------Decoder--------------------------
%%
            %positionIdx = zeros (length(equalizedSignal), 1);
            %equalizedSignal = Vin;
            auxreceivedVoltageSignal = mean(receivedVoltageSignal,2);
            auxreceivedVoltageSignal = repmat(auxreceivedVoltageSignal, 1, NumRx);
            auxreceivedVoltageSignal = abs (auxreceivedVoltageSignal - receivedVoltageSignal);
            
            %auxEqualizedSignal = abs(equalizedSignal);
            
            %Decoding transmission index
            [~ , positionIdx] = max(auxreceivedVoltageSignal,[],2);
            for indChannel = 1:length(positionIdx)
                switch positionIdx (indChannel,1)
                    case 1 %00
                        decodedPosition(indChannel, :) = [0 0];
                    case 2 %01
                        decodedPosition(indChannel, :) = [0 1];
                    case 3 %11
                        decodedPosition(indChannel, :) = [1 1];
                    case 4 %10
                        decodedPosition(indChannel, :) = [1 0];
                end
            end
            %Decoding transmitted symbols
            I = (1 : size(receivedVoltageSignal, 1)) .';
            J = reshape(positionIdx, [], 1);
            linearInd = sub2ind(size(receivedVoltageSignal), I, J);
            y = receivedVoltageSignal(linearInd);
            symbolOutput = y;
%-------------------------Decoder--------------------------

%%

%-------------------------Least-Squares Equalization--------------------------
%%
            unbiasedReceivedVoltageSignal = receivedVoltageSignal - VDC;                   %pra que serve esse sinal?? 
            if LS_EQUALIZATION == 1
                [f] = LSEqualizer(index, pilot, y, numberOfBits, COD_GEOMETRY, LED_POSITION);                              
                symbolOutput = filter(f, 1, y) ; 
            end
%-------------------------Least-Squares Equalization--------------------------
%%
%-------------------------Demodulation--------------------------
%%
            %delay = delta;            
            
            decDemodSymbolAux = pam4HardThreshold(symbolOutput);
            
            [corr,lags] = xcorr(decDemodSymbolAux, y);
            [~,idx] = max(abs(corr));
            delay = abs(lags(idx));
            
            decDemodSymbol = pamdemod(decDemodSymbolAux,2^numberOfBits,0,'gray'); %sinal demodulado               
            decodedSymbol = de2bi(decDemodSymbol,numberOfBits);
%-------------------------Demodulation--------------------------
%%           
            bitOutputData = zeros (length(bitInputData), 2);
            bitOutputData(1:2:end) = [decodedSymbol(delay+1:bitInputDataLength,:) ;zeros(delay,2)];
            bitOutputData(2:2:end) = decodedPosition;
            
            bitOutputDataSymbol = zeros (length(bitInputData), 2);
            bitOutputDataPosition = zeros (length(bitInputData), 2);
            bitOutputDataSymbol(1:2:end) = [decodedSymbol(delay+1:bitInputDataLength,:) ;zeros(delay,2)];
            bitOutputDataPosition(2:2:end) = decodedPosition;
            
            berAuxSymbol(j)     = sum(sum(abs(bitOutputDataSymbol(  1:(blockLength/2),:) - bitInputDataSymbol(1:(blockLength/2),:))))./(blockLength/2);
            berAuxPosition(j)   = sum(sum(abs(bitOutputDataPosition(1:(blockLength/2),:) - bitInputDataPosition(1:(blockLength/2),:))))./(blockLength/2);
            berAux(j)           = sum(sum(abs(bitOutputData(        1:(blockLength/2),:) - bitInputData(1:(blockLength/2),:))))./blockLength;
        end
        berSymbol(SNRIndex,index) = mean(berAuxSymbol);
        berPosition(SNRIndex,index) = mean(berAuxPosition);
        
        ber(SNRIndex,index) = mean(berAux);
    end
end

save(['.' filesep 'resultsBER' filesep FileName '.mat'],'ber','berPosition', 'berSymbol','SNR', 'monteCarloLoops', 'numberOfSymbols', 'modulationIndexVector');