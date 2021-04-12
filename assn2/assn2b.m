clear all;
alpha= 0.35;
A = 5;
T = 1/25000;
Ts = T/8;
t = [-A*T:Ts:A*T] + 10^(-8);
num = sin(pi*(1-alpha)*t/T) + (4*alpha*t/T).*cos(pi*(1+alpha)*t/T);
den = (1-(4*alpha*t/T).^2).*(pi*t/T);
h = num./den;

esN0_dB = [0:2:14];
ps = zeros(1, length(esN0_dB));
pn = zeros(1, length(esN0_dB));

noErrs = zeros(1, length(esN0_dB));
noBitErrs = zeros(1, length(esN0_dB));
ser = zeros(1, length(esN0_dB));
ber = zeros(1, length(esN0_dB));
trails = 100;
for i = 1:length(esN0_dB)
    for j = 1:trails
        QPSK =[-1 1];
        N = 1024;
        x= randsrc(1,N,QPSK) + 1i*randsrc(1,N,QPSK);
        S = (1/sqrt(2))*x;

        y = upsample(S, 8);
        w = conv(y, h);
        ww = conv(w, h);
        
        ps = norm(ww)^2;
        
        noise = 1/sqrt(2)*[randn(1,length(w)) + 1i*randn(1,length(w))];
        noise = 10^(-esN0_dB(i)/20)*noise;
        
        nw = conv(noise, h);
        
        pn = norm(nw)^2;
        %estimate alpha from ps, pn and dB
        alpha1 = sqrt(ps/(pn*(10^((esN0_dB(i)+3)/10))));

        
        final = ww + alpha1 * nw;
        
        rs = final(81:8:end-80);
        rsre = real(rs); % real part
        rsim = imag(rs); % imaginary part

        ip_hat = zeros(1,N);

        ip_hat(find(rsre < 0 & rsim< 0)) = -1 + -1i;
        ip_hat(find(rsre >= 0 & rsim > 0)) = 1 + 1i;
        ip_hat(find(rsre < 0 & rsim >= 0)) = -1 + 1i;
        ip_hat(find(rsre >= 0 & rsim < 0)) = 1 - 1i;

        noErrs(i) =  noErrs(i) + size(find([x- ip_hat]),2);
        noBitErrs(i) = noBitErrs(i) + size(find([real(x)- real(ip_hat)]),2) + size(find([imag(x)- imag(ip_hat)]),2);
        
    end
end

ser = noErrs/(N*trails);
ber = noBitErrs/ (2*N*trails);

figure(2)
semilogy(esN0_dB, ser, 'g')
grid on
ylabel('SER')
xlabel('E_s/N_0 (dB)')
title('Symbol Error Rate for Coherent QPSK')
figure(3)
BER = 1/2.*erfc(sqrt(10.^(esN0_dB/10)));
plot(esN0_dB, ber,'r',esN0_dB, BER, '--')
grid on
hold on
ylabel('BER')
xlabel('E_b/N_0 (dB)')
title('Bit Error Rate for Coherent QPSK')
legend('calculated', 'theoretical')