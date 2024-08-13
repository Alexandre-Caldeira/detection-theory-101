% OBJ: Gerar thresholds CGST-Beta e aplicar em MSC MC de sinal puro e ruido


%% Setup
clearvars; close all; clc

%% Parametros
K = 5;           % numero total de testes a aplicar sequencialmente
FPd = 0.01;      % taxa de falso positivo desejado para o exame
Mmax_exame = 90; % tempo total de exame em segundos
M = 90/K;        % numero janelas de 1s usadas para cada teste

[aThresholds,gThresholds] = get_beta_CGST_thresholds(K, M, FPd);

%% Teste com distribuição teórica

disp('-------------------------------------------------------------------')
disp('MC amostrando Beta')
FP      = zeros(1, K);      % number of false-positives
TN      = zeros(1, K);      % number of true-negatives
NumT    = 1*1e6;          % number of tests to carry out

for ti=1:NumT
    Ps          = betarnd(1,M-1,1,K);            % random p values: all assumptions are met, so the FPR should be exact. 
    Plog        = Ps;           % fisher transform
    Detected    = false;                % 
    
    % check for rejections
    for k=1:K
        if sum(Plog(1:k)) >= aThresholds(k)
            FP(k) = FP(k)+1;
            break
        elseif sum(Plog(1:k)) <= gThresholds(k)
            TN(k) = TN(k)+1;
            break
        end
    end
end

Stage_FPRs  = FP/NumT               % stage-wise FPRs
Stage_TNRs  = TN/NumT               % stage-wise FPRs
FPR         = sum(FP) / NumT        % total FPR
TNR         = sum(TN) / NumT        % total FPR

%% Teste com MSC de ruido
% 
% disp('MC amostrando FFT(rndn)')
% FS = 1000; 
% NumT = 5e3;
% 
% % For each test from 1 to NumT
% % MSC = nan(NumT,FS/2+1);
% % Generate Noise with M windows
% NFFT = FS;                         % rever isso aqui ! motivo/causas/impactos comp. 
% Njanelas = M;
% NpontosTotal = NFFT*Njanelas;      % Numero total de pontos de cada sinal
% ruido = randn(NumT,NpontosTotal);  % Gera um ruido gaussiano, teoricamente de variancia unitaria e media nula
% ruido = ruido-mean(ruido);         % Força a media nula
% signals = ruido./std(ruido);       % Força a variancia desejada para o sinal
% signals = reshape(signals,NumT*NFFT, Njanelas);
% 
% % Compute FFT of all M windows
% SIGNALS = fft(signals);
% SIGNALS = SIGNALS(1:floor(end/2)+1,:); % only half the FFT spectrum is valid
% f = FS/2*linspace(0,1,NFFT/2+1)'; % only half the FFT spectrum is valid
% 
% % Compute MSC in any freq for the M-th window
% MSC = msc_fft(SIGNALS,M);
% MSC = reshape(MSC,[],1);
% 
% % histogram(MSC) % see simulated H0 PDF
% count = 0;
% FP      = zeros(1, K);      % number of false-positives
% TN      = zeros(1, K);      % number of true-negatives
% 
% for ti=1:5:NumT
%     Ps          = MSC(ti:ti+4); 
%     Plog        = Ps;           
%     count = count+1;
%     % check for rejections
%     for k=1:K        
%         if sum(Plog(1:k)) >= aThresholds(k)
%             FP(k) = FP(k)+1;
%             break
%         elseif sum(Plog(1:k)) <= gThresholds(k)
%             TN(k) = TN(k)+1;
%             break
%         end
% 
%     end
% end
% 
% % Stage_FPRs  = FP/count               % stage-wise FPRs
% % Stage_TNRs  = TN/count               % stage-wise FPRs
% FPR         = sum(FP) / count        % total FPR
% TNR         = sum(TN) / count        % total FPR


disp('-------------------------------------------------------------------')
disp('MC amostrando FFT(rndn) SNR = -100 dB')

NumT = 5e3;
FS = 1000;  % frequencia de amostragem 
SFREQ = 80; % frequencia de estimulacao 
NFFT = FS;                         % rever isso aqui ! motivo/causas/impactos comp. 
Nsinais = NumT;
Njanelas = Mmax_exame;

SNRfun = @() -100;
[S5] = gen_signals(SNRfun, FS, SFREQ, NFFT, Nsinais, Njanelas);

% histogram(MSC) % see simulated H0 PDF
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

Stage_FPRs  = FP/count               % stage-wise TPRs
Stage_TNRs  = TN/count               % stage-wise FNRs
FPR         = sum(FP) / count        % total TPR
TNR         = sum(TN) / count        % total FNR
count

%% Teste com MSC de sinal ruidoso (SNR = 15 dB)

disp('-------------------------------------------------------------------')
disp('MC amostrando FFT(sinal+rndn) SNR = 15 dB')
FS = 1000;  % frequencia de amostragem 
SFREQ = 80; % frequencia de estimulacao 
NFFT = FS;                         % rever isso aqui ! motivo/causas/impactos comp. 
Nsinais = NumT;
Njanelas = Mmax_exame;

SNRfun = @() 15;
[S5] = gen_signals(SNRfun, FS, SFREQ, NFFT, Nsinais, Njanelas);

% histogram(MSC) % see simulated H0 PDF
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

Stage_TPRs  = TP/count               % stage-wise TPRs
Stage_FNRs  = FN/count               % stage-wise FNRs
TPR         = sum(TP) / count        % total TPR
FNR         = sum(FN) / count        % total FNR

%% Teste com MSC de sinal ruidoso (SNR = -5 dB)

disp('-------------------------------------------------------------------')
disp('MC amostrando FFT(sinal+rndn) SNR = -5 dB')
FS = 1000;  % frequencia de amostragem 
SFREQ = 80; % frequencia de estimulacao 
NFFT = FS;                         % rever isso aqui ! motivo/causas/impactos comp. 
Nsinais = NumT;
Njanelas = Mmax_exame;

SNRfun = @() -5;
[S5] = gen_signals(SNRfun, FS, SFREQ, NFFT, Nsinais, Njanelas);

% histogram(MSC) % see simulated H0 PDF
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

Stage_TPRs  = TP/count               % stage-wise TPRs
Stage_FNRs  = FN/count               % stage-wise FNRs
TPR         = sum(TP) / count        % total TPR
FNR         = sum(FN) / count        % total FNR

%% Teste com MSC de sinal ruidoso (SNR = -10 dB)

disp('-------------------------------------------------------------------')
disp('MC amostrando FFT(sinal+rndn) SNR = -10 dB')
FS = 1000;  % frequencia de amostragem 
SFREQ = 80; % frequencia de estimulacao 
NFFT = FS;                         % rever isso aqui ! motivo/causas/impactos comp. 
Nsinais = NumT;
Njanelas = Mmax_exame;

SNRfun = @() -10;
[S5] = gen_signals(SNRfun, FS, SFREQ, NFFT, Nsinais, Njanelas);

% histogram(MSC) % see simulated H0 PDF
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

        Plog(k) = Ps(80);

        if sum(Plog(1:k)) >= aThresholds(k)
            TP(k) = TP(k)+1;
            break
        elseif sum(Plog(1:k)) <= gThresholds(k)
            FN(k) = FN(k)+1;
            break
        end

    end
end

% Stage_TPRs  = TP/count               % stage-wise TPRs
% Stage_FNRs  = FN/count               % stage-wise FNRs
TPR         = sum(TP) / count        % total TPR
FNR         = sum(FN) / count        % total FNR

%% Teste com MSC de sinal ruidoso (SNR = -15 dB)

disp('-------------------------------------------------------------------')
disp('MC amostrando FFT(sinal+rndn) SNR = -15 dB')
FS = 1000;  % frequencia de amostragem 
SFREQ = 80; % frequencia de estimulacao 
NFFT = FS;                         % rever isso aqui ! motivo/causas/impactos comp. 
Nsinais = NumT;
Njanelas = Mmax_exame;

SNRfun = @() -15;
[S5] = gen_signals(SNRfun, FS, SFREQ, NFFT, Nsinais, Njanelas);

% histogram(MSC) % see simulated H0 PDF
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

        Plog(k) = Ps(80);

        if sum(Plog(1:k)) >= aThresholds(k)
            TP(k) = TP(k)+1;
            break
        elseif sum(Plog(1:k)) <= gThresholds(k)
            FN(k) = FN(k)+1;
            break
        end

    end
end

% Stage_TPRs  = TP/count               % stage-wise TPRs
% Stage_FNRs  = FN/count               % stage-wise FNRs
TPR         = sum(TP) / count        % total TPR
FNR         = sum(FN) / count        % total FNR

%% Teste com MSC de sinal ruidoso (SNR = -25 dB)

disp('-------------------------------------------------------------------')
disp('MC amostrando FFT(sinal+rndn) SNR = -25 dB')
FS = 1000;  % frequencia de amostragem 
SFREQ = 80; % frequencia de estimulacao 
NFFT = FS;                         % rever isso aqui ! motivo/causas/impactos comp. 
Nsinais = NumT;
Njanelas = Mmax_exame;

SNRfun = @() -25;
[S5] = gen_signals(SNRfun, FS, SFREQ, NFFT, Nsinais, Njanelas);

% histogram(MSC) % see simulated H0 PDF
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

        Plog(k) = Ps(80);

        if sum(Plog(1:k)) >= aThresholds(k)
            TP(k) = TP(k)+1;
            break
        elseif sum(Plog(1:k)) <= gThresholds(k)
            FN(k) = FN(k)+1;
            break
        end

    end
end

% Stage_TPRs  = TP/count               % stage-wise TPRs
% Stage_FNRs  = FN/count               % stage-wise FNRs
TPR         = sum(TP) / count        % total TPR
FNR         = sum(FN) / count        % total FNR






