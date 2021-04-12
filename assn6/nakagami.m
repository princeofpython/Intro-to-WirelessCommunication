snr_db = (0:1:40);
ber = zeros(1, length(snr_db));
m_values = [0.5 ,1, 1.5, 2, 2.5, 3, 4, 5];
figure(1)

for ii = 1:length(m_values)
    ber = (1.+ (10.^(snr_db/10))/m_values(ii)).^(-m_values(ii));
    semilogy(snr_db, ber)
    grid on
    hold on
    ylabel('BER')
    xlabel('E_b/N_0 (dB)')
    title('Bit Error Rate for DBPSK')
end
BER = 1/2.*exp(-10.^(snr_db/10));
semilogy(snr_db, BER)
ylim([10^(-9) 10^0])
legend('m = 0.5', 'm = 1', 'm = 1.5', 'm = 2', 'm = 2.5', 'm = 3','m = 4','m = 5', 'No fading')