function treta = pam4HardThreshold(x)


treta = x;

treta(treta >= 0 & treta<2) = 1;

treta(treta >=2) = 3;

treta(treta<0 & treta>= -2) = -1;

treta(treta < -2) = -3;


