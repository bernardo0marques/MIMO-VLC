function treta = pam8HardThreshold(x)


treta = x;

treta(treta >= 0 & treta<2) = 1;

treta(treta >= 2 & treta<4) = 3;

treta(treta >= 4 & treta<6) = 5;

treta(treta >=6) = 7;

treta(treta<0 & treta>= -2) = -1;

treta(treta<-2 & treta>= -4) = -3;

treta(treta<-4 & treta>= -6) = -5;

treta(treta < -6) = -7;


