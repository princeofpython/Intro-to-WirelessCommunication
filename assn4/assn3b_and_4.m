%This also contains solution for assn 4.
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

esN0_dB = [0:1:10];
ps = zeros(1, length(esN0_dB));
pn = zeros(1, length(esN0_dB));

noBitErrsRaleigh = zeros(1, length(esN0_dB));
noBitErrs = zeros(1, length(esN0_dB));

trails = 2000;
for i = 1:length(esN0_dB)
    for j = 1:trails
        x = randsrc(1, N, DBPSK);
        prev = 1;
        s = zeros(1,N);
        deltatheta = [0, pi];
        for ii = 1:N
            s(ii)= prev*exp(1i*deltatheta(x(ii)+1));
            prev = s(ii);
        end

        y = upsample(s, 8);
        w = conv(y, h);
        W = conv(w, h);
        
        ps = norm(W(81:8:end-80))^2;
        w1 = awgn(w, esN0_dB(i));
        %noise = 10^(-esN0_dB(i)/20)*noise;
        noise = w1 - w;
        nw = conv(noise, h);
        
        pn = norm(nw(81:8:end-80))^2;
        %estimate alpha from ps, pn and dB
        alpha1 = sqrt(ps/(pn*(10^((esN0_dB(i))/10))));

        
        raleighxy = sqrt(0.5).*randn(1,2);
        alphasquareraleigh = sum(raleighxy.^2);
        
        final = W + alpha1 * nw;
        finalraleigh = W + (alpha1/alphasquareraleigh)*nw;
        
        rs = final(81:8:end-80);
        rsrl = finalraleigh(81:8:end-80);
        rx = zeros(1,N);
        rxrl = zeros(1,N);
        prev = 1;
        prevrl = 1;
        for ii = 1:N
            if real(rs(ii)*conj(prev))>0
                rx(ii)= 0;
            else
                rx(ii) = 1;
            end
            if real(rsrl(ii)*conj(prevrl))>0
                rxrl(ii)= 0;
            else
                rxrl(ii) = 1;
            end

            prev = rs(ii);
            prevrl = rsrl(ii);
        end
        noBitErrsRaleigh(i) = noBitErrsRaleigh(i) + size(find([x- rxrl]),2);
        noBitErrs(i) = noBitErrs(i) + size(find([x- rx]),2);
    end
    

end
raleigherrors = zeros(1, length(esN0_dB));
for i = 1:length(esN0_dB)
    for j = 1:1000
        raleighxy = sqrt(0.5).*randn(1,2);
        alphasquareraleigh = sum(raleighxy.^2);
        
        raleigherrors(i) = raleigherrors(i)+1/2*exp(-alphasquareraleigh*10^(esN0_dB(i)/10));
        
    end
end
berraleigh = raleigherrors/1000;
berrl = noBitErrsRaleigh/(N*trails);
ber = noBitErrs/ (N*trails);
figure(1)
BER = 1/2.*exp(-10.^(esN0_dB/10));
semilogy(esN0_dB, ber,'r')

grid on
hold on
semilogy(esN0_dB, BER,'b')
semilogy(esN0_dB, berraleigh,'g')
BERRL = 1./(2*(1+(10.^(esN0_dB/10))));
semilogy(esN0_dB, BERRL,'k')
ylabel('BER')
xlabel('E_b/N_0 (dB)')
title('Bit Error Rate for DBPSK')
legend('calculated', 'theoretical', 'RF calculated','RF theoretical')