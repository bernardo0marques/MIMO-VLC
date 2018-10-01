% Comparação dos IM - Equalizado x Não equalizado
addpath(['.' filesep 'resultsBER']);
addpath(['.' filesep 'plots']);
figProp = struct('size',38,'font','Times','lineWidth',2,'figDim',[1 1 600 400]);

for CONSTELATION_SIZE = 2:2
    for LED_POSITION = 1:1
        for TECHNIQUE = 3:3
            for LS_EQUALIZATION = 0:0 
                LED_GEOMETRY = 0; 
                if (LED_POSITION > 2), LED_GEOMETRY = 1; end
                CONFIG;
                load ([FileName '.mat']);
                figure;
                for idx = 1:length(modulationIndexVector)
                    semilogy(SNR, ber(:,idx));
                    hold on;
                end
                x1 = pi;
                y1 = sin(pi);
                txt1 = ' 2\times10^{-1}';
                text(20,0.3e-1,txt1)

                legend({['IM = $', num2str(modulationIndexVector(1)), '$']...
                       ,['IM = $', num2str(modulationIndexVector(2)), '$']...
                       ,['IM = $', num2str(modulationIndexVector(3)), '$'],},  'Location','southwest','Interpreter', 'LaTex')
                xlim([0,30]);
                ylim([1E-6,1]);
                xlabel('SNR [dB]');
                ylabel('BER');
                hold off;
                l.FontSize = 38;
                figFileName = ['.' filesep 'plots' filesep FileName];
                formatFig(gcf,figFileName,'pt',figProp);
            end
        end
    end
end
close all;


