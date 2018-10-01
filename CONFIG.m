%-------------------------String Codes--------------------------
if CONSTELATION_SIZE == 1, M = 4; end
if CONSTELATION_SIZE == 2, M = 8; end

if M == 4, COD_SYMBOL = '4PAM'; end
if M == 8, COD_SYMBOL = '8PAM'; end

if TECHNIQUE == 1, COD_TECHNIQUE = 'RC'; end
if TECHNIQUE == 2, COD_TECHNIQUE = 'SMP'; end
if TECHNIQUE == 3, COD_TECHNIQUE = 'SM'; end

if LS_EQUALIZATION == 0, COD_EQUALIZATION = 'NEQ'; end
if LS_EQUALIZATION == 1, COD_EQUALIZATION = 'EQ'; end

if LED_GEOMETRY == 0, COD_GEOMETRY = 'CIRC'; end

if LED_POSITION == 1, R_tx = 10e-2; end 
if LED_POSITION == 1, R_rx = 10e-2; end 
if LED_POSITION == 1, d    = 20e-2; end 

if LED_POSITION == 2, R_tx = 10e-2; end 
if LED_POSITION == 2, R_rx = 10e-2; end 
if LED_POSITION == 2, d    = 175e-2; end 

if LED_POSITION == 3, R_tx = 3e-2; end 
if LED_POSITION == 3, R_rx = R_tx; end
if LED_POSITION == 3, d    = 20e-2; end

if LED_GEOMETRY == 1, COD_GEOMETRY = 'LIN'; end

FileName = [COD_SYMBOL, '_', COD_TECHNIQUE, '_', COD_EQUALIZATION, '_', COD_GEOMETRY, '_P', int2str(LED_POSITION)];
ScriptName = ['VLC_',COD_SYMBOL, '_', COD_TECHNIQUE,'.m'];
%-------------------------String Codes--------------------------