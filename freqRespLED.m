function G = freqRespLED(w)

wc = 2*pi*1e6; % 1,0 MHz x 2pi = 6,283 Mrad/s

w2 = 2*pi*3.26e6; %  2,6MHz x 2pi = 20,48 Mrad/s

w3 = 2*pi*10.86e6; % 8,6MHz x 2pi = 68,24 Mrad/s


wAux1 = w(abs(w)<=wc); 

wAux2 = w(abs(w)>wc);

G = zeros(length(wAux1) + length(wAux2),1);

G(abs(w)<=wc) = exp(-abs(wAux1)/w2);

G(abs(w)>wc) = exp(-abs(wc)/w2)*exp(abs(wc)/w3)*exp(-abs(wAux2)/w3);

