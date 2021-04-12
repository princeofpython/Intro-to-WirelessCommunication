clear all;
alpha= 0.35;
    A = 5;
    T = 1/25000;
    Ts = T/8;
    t = [-A*T:Ts:A*T] + 10^(-8);
    num = sin(pi*(1-alpha)*t/T) + (4*alpha*t/T).*cos(pi*(1+alpha)*t/T);
    den = (1-(4*alpha*t/T).^2).*(pi*t/T);
    h1 = num./den;
    figure(1)
    plot(t/T,h1,'LineWidth',1.5);
    hold on
alpha= 0.7;
    num = sin(pi*(1-alpha)*t/T) + (4*alpha*t/T).*cos(pi*(1+alpha)*t/T);
    den = (1-(4*alpha*t/T).^2).*(pi*t/T);
    h2 = num./den;
    plot(t/T,h2,'LineWidth',1.5);
    hold on
alpha= 1;
    num = sin(pi*(1-alpha)*t/T) + (4*alpha*t/T).*cos(pi*(1+alpha)*t/T);
    den = (1-(4*alpha*t/T).^2).*(pi*t/T);
    h3 = num./den;
    plot(t/T,h3,'LineWidth',1.5);
    hold on
grid on
xlabel('t/T (Seconds)')
ylabel('Amplitude')
title('Time Domain SRRC pulse')
xticks([-5 -4 -3 -2 -1 0 1 2 3 4 5])
legend('\alpha = 0.35','\alpha = 0.7','\alpha = 1')
xlim([-5 5])

figure(2);

H1 = fft(h1,500);
H_norm1 = abs(H1)/max(abs(H1));
H2 = fft(h2,500);
H_norm2 = abs(H2)/max(abs(H2));
H3 = fft(h3,500);
H_norm3 = abs(H3)/max(abs(H3));
f = 1:length(H_norm1)/2;
plot(f/max(f),20*log10(H_norm1(1:length(H_norm1)/2)),'LineWidth',1.5);
hold on
plot(f/max(f),20*log10(H_norm2(1:length(H_norm2)/2)),'LineWidth',1.5);
plot(f/max(f),20*log10(H_norm3(1:length(H_norm3)/2)),'LineWidth',1.5);
grid on;
title('Magnitude Response of SRRC pulse')
legend('\alpha = 0.35','\alpha = 0.7','\alpha = 1')
xlabel('f in Hz')

w1 = conv(h1,h1);
w2 = conv(h2,h2);
w3 = conv(h3,h3);
tf = [-2*A*T:Ts:2*A*T];
figure(3);
plot(tf,w1,'LineWidth',1.5);
hold on;
plot(tf,w2,'LineWidth',1.5);
plot(tf,w3,'LineWidth',1.5);
grid on;
legend('\alpha = 0.35','\alpha = 0.7','\alpha = 1')
xlabel('t/T (Seconds)')
ylabel('Amplitude')
title('Time Domain RC pulse')

signls = [ -1 -1 1 1 1 -1 1 1 -1 -1 -1 1 -1 1 -1 1 -1 1 1 1];

taxis = [-10*T:Ts:32*T];
finalwaveform1 = zeros(1,length(taxis));
finalwaveform2 = zeros(1,length(taxis));
finalwaveform3 = zeros(1,length(taxis));
for i=1:length(signls)
    finalwaveform1(9*i-8:152+9*i) = finalwaveform1(9*i-8:152+9*i) + signls(i)* w1;
    finalwaveform2(9*i-8:152+9*i) = finalwaveform2(9*i-8:152+9*i) + signls(i)* w2;
    finalwaveform3(9*i-8:152+9*i) = finalwaveform3(9*i-8:152+9*i) + signls(i)* w3;
end
figure(4)
plot(taxis/T, finalwaveform1,'LineWidth',1.5);
xlabel('t/T (Seconds)')
ylabel('Amplitude')
title('Pulse shaped pulse for \alpha = 0.35')
grid on;
figure(5)
plot(taxis/T, finalwaveform2,'LineWidth',1.5);
xlabel('t/T (Seconds)')
ylabel('Amplitude')
title('Pulse shaped pulse for \alpha = 0.7');
grid on;
figure(6)
plot(taxis/T, finalwaveform3,'LineWidth',1.5);
xlabel('t/T (Seconds)')
ylabel('Amplitude')
title('Pulse shaped pulse for \alpha = 1.0')
grid on;

rxsignls= [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
symbols = [ 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
for i=1:20
    symbols(i) = finalwaveform1(72+8*i);
    if finalwaveform1(72+8*i)>=0
        rxsignls(i) = 1;
    else
        rxsignls(i) = -1;
    end
end
