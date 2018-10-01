tic

clear;clc;close all;
addpath(['..' filesep 'Degree_Radian_Conversion']);
addpath(['.' filesep 'LED Parameters']);
addpath(['.' filesep 'resultsBER']);
addpath(['.' filesep 'Equalizer']);

load whiteLED_334-15.mat;

%%
%-------------------------Simulation Scenarios--------------------------
monteCarloLoops = 10;
numberOfSymbols = 1000;
%-------------------------Simulation Scenarios--------------------------
%TECHNIQUE          (1/2/3)   - RC/SMP/SM
%LS_EQUALIZATION    (0/1)     - Não Equalizado/Equalizado
%LED_GEOMETRY       (0/1)     - Circular/Em linha 
%LED_POSITION       (1/2/3)   - R_rx = ((20cm/40cm/60cm)/R_tx = 3cm) 
%-------------------------Simulation Scenarios--------------------------
%%
for CONSTELATION_SIZE = 1:2
    for LED_POSITION = 1:3
        for TECHNIQUE = 3:3
            for LS_EQUALIZATION = 1:1 
                LED_GEOMETRY = 0; 
                if (LED_POSITION > 2), LED_GEOMETRY = 1; end
                CONFIG;  
                VLC_Define;      
                run(ScriptName);
            end
        end
    end
end

%compare_BER_IM;
%compare_BER_Techniques;

toc
minutos = toc/60 
load chirp, sound(y,1/2*Fs)
