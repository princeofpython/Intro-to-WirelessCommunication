clear all;
alpha= 0.35;
A = 5;
T = 1/25000;
Ts = T/8;
t = [-A*T:Ts:A*T] + 10^(-8);
num = sin(pi*(1-alpha)*t/T) + (4*alpha*t/T).*cos(pi*(1+alpha)*t/T);
den = (1-(4*alpha*t/T).^2).*(pi*t/T);
h = num./den;

DBPSK = [0 1]; 
N = 512;
x = randsrc(1, N, DBPSK);
prev = 1;
s = zeros(1,N);
deltatheta = [0, pi];
for i = 1:N
    s(i)= prev*exp(1i*deltatheta(x(i)+1));
    prev = s(i);
end

y = upsample(s, 8);
w = conv(y, h);

w1 = awgn(w, 6-10*log10(8));
W = conv(w1, h);
rs = W(81:8:end-80);
rsre = real(rs); % real part
rsim = imag(rs); % imaginary part
figure(1)
scatter(real(rs), imag(rs));


rx = zeros(1,N);
prev = 1;
for i = 1:N
    if real(rs(i)*conj(prev))>0
        rx(i)= 0;
    else
        rx(i) = 1;
    end
        
    prev = rs(i);
end
noErr = size(find([x- rx]),2);



