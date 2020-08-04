function [ R ] = coherence( v, v2, Fs )
    % Constants
    PLI_S = 56;      % Power line interference frequency                          
    PLI_E = 64;                                                                   
    PL2_S = (Fs-180)-2;     % Power line interference frequency second band       
    PL2_E = (Fs-180)+2;                                                           
    PL3_S = 117;      % Power line interference frequency third band              
    PL3_E = 123; 
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

    [cxx,f] = mscohere(v', v2',hamming(2*Fs),[],[],Fs);
    cxx = sqrt(cxx);

    mask_del = ((f > DEL_S) & (f < DEL_E));
    mask_the = ((f >= THE_S) & (f < THE_E));
    mask_alp = ((f >= ALP_S) & (f < ALP_E));
    mask_bet = ((f >= BET_S) & (f < THZ_S)) | ((f > THZ_E) & (f < BET_E));
    %mask_bet = ((f >= BET_S) & (f < BET_E));
    mask_gam = ((f >= GAM_S) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < GAM_E));
    %mask_bro = ((f > BRO_S) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < PL3_S)) | ((f > PL3_E) & (f < BRO_E));
    mask_bro = ((f > BRO_S) & (f < THZ_S)) | ((f > THZ_E) & (f < PLI_S)) | ((f > PLI_E) & (f < PL2_S)) | ((f > PL2_E) & (f < PL3_S)) | ((f > PL3_E) & (f < BRO_E));
    cohDel = mean(cxx(mask_del));
    cohThe = mean(cxx(mask_the));
    cohAlp = mean(cxx(mask_alp));
    cohBet = mean(cxx(mask_bet));
    cohGam = mean(cxx(mask_gam));
    cohBro = mean(cxx(mask_bro));
    
    R = zeros(6,1);
    R(1) = cohDel;
    R(2) = cohThe;
    R(3) = cohAlp;
    R(4) = cohBet;
    R(5) = cohGam;
    R(6) = cohBro;

end

