clear all;
alpha= 0.35;
A = 5;
T = 1/25000;
Ts = T/8;
t = [-A*T:Ts:A*T] + 10^(-8);
num = sin(pi*(1-alpha)*t/T) + (4*alpha*t/T).*cos(pi*(1+alpha)*t/T);
den = (1-(4*alpha*t/T).^2).*(pi*t/T);
h = num./den;

QPSK = [-1 1]; 
N = 512;
x= randsrc(1,N,QPSK) + 1i*randsrc(1,N,QPSK);
S = (1/sqrt(2))*x;

y = upsample(S, 8);
w = conv(y, h);



noise = 1/sqrt(2)*[randn(1,length(w)) + 1i*randn(1,length(w))];
noise = 10^(-6/20)*noise;

w1 = w+ noise;

W = conv(w1, h);

rs = W(81:8:end-80);
rsre = real(rs); % real part
rsim = imag(rs); % imaginary part

ip_hat = zeros(1,N);

ip_hat(find(rsre < 0 & rsim< 0)) = -1 + -1*j;
ip_hat(find(rsre >= 0 & rsim > 0)) = 1 + 1*j;
ip_hat(find(rsre < 0 & rsim >= 0)) = -1 + 1*j;
ip_hat(find(rsre >= 0 & rsim < 0)) = 1 - 1*j;



noErr = size(find([x- ip_hat]),2);
figure(1)
scatter(real(rs), imag(rs));


