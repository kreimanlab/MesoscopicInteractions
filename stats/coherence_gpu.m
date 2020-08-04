function [ R ] = coherence_gpu( v, nperseg, fs, pli )

    % Inputs
    %   v - matrix of timeseries organized by columns
    % Constants
    PLI_S = pli-4;      % Power line interference frequency                          
    PLI_E = pli+4;                                                                   
    PL2_S = (fs-(3*pli))-2;     % Power line interference frequency second band       
    PL2_E = (fs-(3*pli))+2;                                                           
    PL3_S = (2*pli)-3;      % Power line interference frequency third band              
    PL3_E = (2*pli)+3; 
    %
    THZ_S = 17;
    THZ_E = 23;
    %
    DEL_S = 0.5;     % Delta wave                                                 
    DEL_E = 3;                                                                    
    THE_S = 3;       % Theta wave                                                 
    THE_E = 8;                                                                    
    ALP_S = 8;       % Alpha wave                                                 
    ALP_E = 12;                                                                   
    BET_S = 12;      % Beta wave                                                  
    BET_E = 30;                                                                   
    GAM_S = 30;      % Gamma wave                                                 
    GAM_E = 100;                                                                  
    BRO_S = 0.5;     % Broadband                                                  
    BRO_E = 125;
    N_BANDS = 6;
    %--------------------------------------------------

    %[cxx,f] = mscohere(v', v2',hamming(round(fs/5)),round(0.8*round(fs/5)),[],round(fs/5));
    %cxx = sqrt(cxx);
    
    [n_samples,n_chan] = size(v);
    
    % Segmentation
    noverlap = round(nperseg * 0.8);
    stride = nperseg - noverlap;
    n_seg = floor((n_samples - noverlap) / stride);

    %  -=- Zero-pad end of segment MATLAB-style -=-
    nfft = round(max([256,2^(ceil(log2(nperseg)))]));
    % -=- Don't zero-pad Scipy-style -=-
    nfft = nperseg;
    
    % Init segmented matrix
    V = zeros(nfft,n_chan,n_seg);
    f = linspace(0,(fs/2),ceil(nfft/2)+1);
    
    
    mask_del = ((f > DEL_S) & (f < DEL_E));
    mask_the = ((f >= THE_S) & (f < THE_E));
    mask_alp = ((f >= ALP_S) & (f < ALP_E));
    mask_bet = ((f >= BET_S) & (f < THZ_S)) | ((f > THZ_E) & (f < BET_E));
    %mask_bet = ((f >= BET_S) & (f < BET_E));
    mask_gam = ((f >= GAM_S) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < GAM_E));
    %mask_bro = ((f > BRO_S) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < PL3_S)) | ((f > PL3_E) & (f < BRO_E));
    mask_bro = ((f > BRO_S) & (f < THZ_S)) | ((f > THZ_E) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < PL3_S)) | ((f > PL3_E) & (f < BRO_E));
    
    % -=- Hann window (Scipy default) -=-
    %win = hanning(nperseg);
    %  -=- Hamming window (MATLAB default), more robust to PLI -=-
    win = hamming(nperseg);
    
    % Used to normalize out the effect of the window in the final power
    % (this does not affect the coherence)
    %win_s = numpy.sum(win ** 2)
    
    % Loop through segments
    for i = 1:n_seg
        start_i = (i-1)*stride + 1;
        end_i = start_i + nperseg - 1;
        % pad = 0;
        % Apply window
        for j = 1:n_chan
            V(1:nperseg,j,i) = (v(start_i:end_i,j)).*win;
        end
    end
    
    % FFT
    % gpuArray(V)
    Vf = fft(V);
    
    % Get reals
    Vf = Vf(1:(ceil(nfft/2)+1),:,:);
    
    
    
    % Coherence
    n_comb = round(0.5*n_chan*(n_chan-1));
    r = zeros(N_BANDS,n_comb);
    rf = zeros(n_comb,ceil(nfft/2)+1);
    
    % Raw CSD
    cf = zeros(n_comb,ceil(nfft/2)+1);
    
    % Coherence phase
    ph = zeros(N_BANDS,n_comb);
    phf = zeros(n_comb,ceil(nfft/2)+1);
    
    Vf_ms = sqrt(sum((abs(Vf)).^2,3));
    Vf_c = conj(Vf);
    c_comb = 1;
    for i = 1:(n_chan-1)
        for j = (i+1):n_chan
            % CSD
            Xii = Vf_ms(:,(i));
            Xjj = Vf_ms(:,(j));
            Xij = sum(squeeze(Vf(:,i,:)) .* squeeze(Vf_c(:,j,:)),2);
            CSD = Xij./(Xii.*Xjj);
            
            % magnitude
            C2 = abs(CSD);
            % phase
            PHI = atan2(imag(CSD),real(CSD));
            
            % save
            cf(c_comb,:) = CSD;
            rf(c_comb,:) = C2;
            phf(c_comb,:) = PHI;
            
            c_comb = c_comb + 1;
            %return
        end
    end
    
    
    %return
    
%     cohDel = abs(mean(cf(:,mask_del),2));
%     cohThe = abs(mean(cf(:,mask_the),2));
%     cohAlp = abs(mean(cf(:,mask_alp),2));
%     cohBet = abs(mean(cf(:,mask_bet),2));
%     cohGam = abs(mean(cf(:,mask_gam),2));
%     cohBro = abs(mean(cf(:,mask_bro),2));
    
    cohDel = mean(abs(cf(:,mask_del)),2);
    cohThe = mean(abs(cf(:,mask_the)),2);
    cohAlp = mean(abs(cf(:,mask_alp)),2);
    cohBet = mean(abs(cf(:,mask_bet)),2);
    cohGam = mean(abs(cf(:,mask_gam)),2);
    cohBro = mean(abs(cf(:,mask_bro)),2);
    
    
    %disp(size(mean(cf(:,mask_the),2)))
    
    
    R = zeros(N_BANDS,n_comb);
    R(1,:) = (cohDel');
    R(2,:) = (cohThe');
    R(3,:) = (cohAlp');
    R(4,:) = (cohBet');
    R(5,:) = (cohGam');
    R(6,:) = (cohBro');

end