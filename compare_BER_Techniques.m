% Comparação da BER entre as 3 técnicas, para a cada índice de modulação 
addpath(['.' filesep 'resultsBER']);addpath(['.' filesep 'plots']);
figProp = struct('size',38,'font','Times','lineWidth',2,'figDim',[1 1 600 400]);

for CONSTELATION_SIZE = 1:2
    for LED_POSITION = 1:3
        for LS_EQUALIZATION = 0:1
            LED_GEOMETRY = 0; 
            if (LED_POSITION > 2), LED_GEOMETRY = 1; end
            for idx = 1:length(modulationIndexVector)
                figure;
                for TECHNIQUE = 1:3
                    CONFIG;
                    load ([FileName '.mat']);
                    semilogy(SNR, ber(:,idx));
                    hold on;
                end
                xlim([0,30]);
                ylim([1E-6,1]);
                xlabel('SNR [dB]')
                ylabel('BER')
                legend('RC','SMP','SM','Location','southwest')
                hold off;

                l.FontSize = 32;
                figFileNameAux = strrep([filesep 'plots' filesep COD_SYMBOL, '_', COD_EQUALIZATION, '_MI', num2str(modulationIndexVector(idx)), '_', COD_GEOMETRY, '_P', int2str(LED_POSITION)],'.','_');
                figFileName = ['.' figFileNameAux];
                formatFig(gcf,figFileName,'pt',figProp);
            end
        end
    end
end
close all;


