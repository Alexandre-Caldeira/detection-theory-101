% OBJ: Gerar thresholds CGST-Beta e aplicando para vizualizar PD-SNR 

%% Setup
clearvars; close all; clc

%% Parametros
K = 5;           % numero total de testes a aplicar sequencialmente
FPd = 0.01;      % taxa de falso positivo desejado para o exame
Mmax_exame = 90; % tempo total de exame em segundos
M = 90/K;        % numero janelas de 1s usadas para cada teste

[aThresholds,gThresholds] = get_beta_CGST_thresholds(K, M, FPd);

NumT = 5e3;
FS = 1000;  % frequencia de amostragem 
SFREQ = 80; % frequencia de estimulacao 
NFFT = FS;                         % rever isso aqui ! motivo/causas/impactos comp. 
Nsinais = NumT;
Njanelas = Mmax_exame;
% snr_vec = 0:-1:-30;
snr_vec = 0:-1:-30;

Stage_FPRs  = zeros(numel(snr_vec),K);        % stage-wise TPRs
Stage_TNRs  = zeros(numel(snr_vec),K);        % stage-wise FNRs
FPR         = zeros(1,numel(snr_vec));        % total TPR
TNR         = zeros(1,numel(snr_vec));        % total FNR

Stage_TPRs  = zeros(numel(snr_vec),K);        % stage-wise TPRs
Stage_FNRs  =  zeros(numel(snr_vec),K);       % stage-wise FNRs
TPR         = zeros(1,numel(snr_vec));        % total TPR
FNR         = zeros(1,numel(snr_vec));        % total FNR

%% Teste com MSC de sinal ruidoso (SNR = 0 a -30 dB)

disp('-------------------------------------------------------------------')
disp('MC amostrando FFT(rndn) SNR = 0 a -30 dB')
tempo = tic();
for idx = 1:numel(snr_vec)
    SNRfun = @() snr_vec(idx);
    [S5] = gen_signals(SNRfun, FS, SFREQ, NFFT, Nsinais, Njanelas);
    
    % get TP, FN
    count = 0;
    TP      = zeros(1, K);      % number of false-positives
    FN      = zeros(1, K);      % number of true-negatives
    
    for ti=1:NumT
        % Ps          = msc_fft(S5(:,:,ti),M);
        % Plog        = Ps;           
        count = count+1;
        % check for rejections
        Plog = zeros(1,K);
    
        for k=1:K        
    
            ind_inicial = M*(k-1)+1;
            ind_final = ind_inicial+M-1;
            Ps = msc_fft(S5(:,ind_inicial:ind_final,ti),M);
    
            Plog(k) = Ps(SFREQ);
    
            if sum(Plog(1:k)) >= aThresholds(k)
                TP(k) = TP(k)+1;
                break
    
            elseif sum(Plog(1:k)) <= gThresholds(k)
                FN(k) = FN(k)+1;
                break
    
            end
    
        end
    end
    
    Stage_TPRs(idx,:)  = TP/count;               % stage-wise TPRs
    Stage_FNRs(idx,:)  = FN/count;               % stage-wise FNRs
    TPR(idx)         = sum(TP) / count;        % total TPR
    FNR (idx)        = sum(FN) / count;        % total FNR

    % get FP, TN
    count = 0;
    FP      = zeros(1, K);      % number of false-positives
    TN      = zeros(1, K);      % number of true-negatives
    
    for ti=1:NumT
        % Ps          = msc_fft(S5(:,:,ti),M);
        % Plog        = Ps;           
        count = count+1;
        % check for rejections
        Plog = zeros(1,K);
    
        for k=1:K        
    
            ind_inicial = M*(k-1)+1;
            ind_final = ind_inicial+M-1;
            Ps = msc_fft(S5(:,ind_inicial:ind_final,ti),M);
    
            Plog(k) = Ps(140);
    
            if sum(Plog(1:k)) >= aThresholds(k)
                FP(k) = FP(k)+1;
                break
    
            elseif sum(Plog(1:k)) <= gThresholds(k)
                TN(k) = TN(k)+1;
                break
    
            end
    
        end
    end
    
    Stage_FPRs(idx,:)   = FP/count;               % stage-wise TPRs
    Stage_TNRs(idx,:)   = TN/count;               % stage-wise FNRs
    FPR(idx)         = sum(FP) / count;        % total TPR
    TNR(idx)         = sum(TN) / count;        % total FNR

    if rem(snr_vec(idx),5)==0
        Finalizou = idx
        tempo_gasto = toc(tempo)
    end
end

% tempo_gasto =
% 
%   647.4677


%% Viz
figure
subplot(2,2,1)
plot(snr_vec, 100*TPR, 'LineWidth',1.2)
hold on
plot(snr_vec, 100*TPR, 'k.','MarkerSize',8)
grid on
ylim([0,100])
title('TPR')
ylabel('PD [%]', 'FontSize',15)
xlabel('SNR [dB]', 'FontSize',15)

subplot(2,2,2)
plot(snr_vec, 100*FPR, 'LineWidth',1.2)
hold on
plot(snr_vec, 100*FPR, 'k.','MarkerSize',8)
grid on
ylim([0,2])
title('FPR')
ylabel('\alpha [%]', 'FontSize',15)
xlabel('SNR [dB]', 'FontSize',15)

subplot(2,2,3)
plot(snr_vec, 100*TNR, 'LineWidth',1.2)
hold on
plot(snr_vec, 100*TNR, 'k.','MarkerSize',8)
grid on
ylim([98,100])
title('TNR')
ylabel('\beta [%]', 'FontSize',15)
xlabel('SNR [dB]', 'FontSize',15)

subplot(2,2,4)
plot(snr_vec, 100*FNR, 'LineWidth',1.2)
hold on
plot(snr_vec, 100*FNR, 'k.','MarkerSize',8)
grid on
ylim([0,100])
title('FNR')
ylabel('Type-II Error [%]', 'FontSize',15)
xlabel('SNR [dB]', 'FontSize',15)
