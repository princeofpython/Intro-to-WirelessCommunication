clear all;
alpha= 0.35;
A = 5;
T = 1/25000;
Ts = T/8;
t = [-A*T:Ts:A*T] + 10^(-8);
num = sin(pi*(1-alpha)*t/T) + (4*alpha*t/T).*cos(pi*(1+alpha)*t/T);
den = (1-(4*alpha*t/T).^2).*(pi*t/T);
h = num./den;

DQPSK = [0 1 2 3]; 
N = 512;
x = randsrc(1, N, DQPSK);
bits = zeros(1, 2*N);
for i = 1:N
    if x(i)==0
        bits(2*i-1)= 0;
        bits(2*i) = 0;
    elseif x(i)==1
        bits(2*i-1)= 0;
        bits(2*i) = 1;
    elseif x(i)==2
        bits(2*i-1)= 1;
        bits(2*i) = 0;
    else
        bits(2*i-1)= 1;
        bits(2*i) = 1;     
    end
end
prev = exp(1i*pi/4);
s = zeros(1,N);
deltatheta = [pi/4, 3*pi/4, -pi/4, -3*pi/4];
for i = 1:N
    s(i)= prev*exp(1i*deltatheta(x(i)+1));
    prev = s(i);
end

y = upsample(s, 8);
w = conv(y, h);

w1 = awgn(w, 6+10*log10(2)-10*log10(8));
W = conv(w1, h);
rs = W(81:8:end-80);
rsre = real(rs); % real part
rsim = imag(rs); % imaginary part
figure(1)
scatter(real(rs), imag(rs));

rx = zeros(1,2*N);
prev = exp(1i*pi/4);
for i = 1:N
    if real(rs(i)*conj(prev))>0
        rx(2*i)= 0;
    else
        rx(2*i) = 1;
    end
    if imag(rs(i)*conj(prev))>0
        rx(2*i-1)= 0;
    else
        rx(2*i-1) = 1;
    end    
    prev = rs(i);
end
noErr = size(find([bits- rx]),2);