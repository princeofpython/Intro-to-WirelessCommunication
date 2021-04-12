function Z=Jakes(fd, Tsample, t0, Num)
% Z=Jakes() - returns fading waveform generated using Jakes’ model
% fd= maximum Doppler frequency
% Tsample= Time resolution of signal in seconds
% t0 = starting time of oscillators
% Num=Number of samples of fading to be generated
N0=20; %Number of oscillators
N=4*N0+2; %N is even but not multiple of 4
scale1=1/sqrt(2*N0); %Scale factor for real part of fading channel
scale2=1/sqrt(2*(N0+1)); %Scale factor for imag part of fading channel
n=1:N0; %Oscillator index
t = (0: Tsample:(Num-1) *Tsample) + t0; %Time duration
fn = fd*cos(2*pi*n/N); %fn – frequencies of oscillators in Jakes’ model
alp=0;
beta=n*pi/(N0+1);
phi=2*pi*rand(N0,1); %Random phases for each oscillator
ZC1=zeros(1, length(t));
ZS1=zeros(1, length(t));
for k=1:N0
 %Real part
 ZC1=ZC1+2*cos(beta(k))*cos(2*pi*fn(k)*t+phi(k));
 %Imag part
 ZS1=ZS1+2*sin(beta(k))*cos(2*pi*fn(k)*t+phi(k));
end
ZC=(ZC1+sqrt(2)*cos(alp)*cos(2*pi*fd*t))*scale2;
ZS=(ZS1+sqrt(2)*sin(alp)*cos(2*pi*fd*t))*scale1;

Z=complex(ZC,ZS);
%figure(1)
%plot(t, abs(Z));
end



