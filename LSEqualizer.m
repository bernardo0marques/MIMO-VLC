function [LScoeff] = LSEqualizer(modulationIndex, pilot, inputSignal, numberOfBits, COD_GEOMETRY, LED_POSITION, varargin) %
    
if numberOfBits == 2, COD_SYMBOL = '4PAM'; end
if numberOfBits == 3, COD_SYMBOL = '8PAM'; end

fileName = strcat(COD_SYMBOL,'_MI', int2str(modulationIndex), '_P', int2str(LED_POSITION), '_Len', int2str(length(inputSignal)), '.mat');
if exist(['.' filesep 'Equalizer' filesep fileName], 'file')
    load (['.' filesep 'Equalizer' filesep fileName]);
else
    
    equalizerLengthLIMIT = 50;
    plotEqualizerLength = zeros (1, equalizerLengthLIMIT);

    for equalizerLength = 1:equalizerLengthLIMIT
        n = equalizerLength - 1;

        Jmin = zeros(1, n+1);
        LScoeffArray = zeros(equalizerLength, n+1);
        for delta = 0:n
            pTraining = length(pilot) - delta ;
            matrizR = toeplitz(inputSignal(n+1:pTraining), inputSignal(n+1:-1:1));  % build matrix R
            Saux = pilot(1:length(pilot));
            S = Saux(n + 1 - delta :pTraining - delta);  

            LScoeffArray(:, delta+1) =  (matrizR' *matrizR)\matrizR' * S;
            Jmin(:, delta + 1) = S'*S-S'*matrizR/ (matrizR'*matrizR)*matrizR' * S;                %Jmin = S'*S-S'*matrizR* inv(matrizR'*matrizR)*matrizR' * S; % Jmin for this f and delta % equalizer is a filter
        end
        plotEqualizerLength(equalizerLength) = min(Jmin);
        %meanJmin(equalizerLength, 1) = sum(Jmin)/size(Jmin,2); 
    end
    %plot(plotEqualizerLength);    
    %ylim([0,30]);
    
    equalizerLength = 5;
    n = equalizerLength - 1;
    delta = round(equalizerLength/2);
    
    pTraining = length(pilot) - delta ;
    matrizR = toeplitz(inputSignal(n+1:pTraining), inputSignal(n+1:-1:1));  % build matrix R
    Saux = pilot(1:length(pilot));
    S = Saux(n + 1 - delta :pTraining - delta);  

    LScoeff =  (matrizR' *matrizR)\matrizR' * S;
       
    EqualizerLenFileName = ['.' filesep 'Equalizer' filesep, 'EqualizerLen','_', COD_SYMBOL,'_P', int2str(LED_POSITION), '_MI', int2str(modulationIndex), '.mat'];
    save(EqualizerLenFileName, 'plotEqualizerLength');    
    save(['.' filesep 'Equalizer' filesep fileName], 'LScoeff', 'delta');
end

