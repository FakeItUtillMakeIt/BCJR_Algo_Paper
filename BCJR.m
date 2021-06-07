function BCJR
% 输入，信息位和校验位
t = [1 1 -1 1 1 -1 1 1 -1 1 1 1 -1 1 -1 -1];
U = [+0.213  -0.371   +0.139  +0.514  +0.539  -0.422  +1.533  +1.457  +0.323   +2.028 -0.414  +1.482   -1.701  -0.175  -0.862  -0.918];
C = [+0.364   0.000   +1.818   0.000  +0.388   0.000  +0.267   0.000  +1.103   0.000  +3.560   0.000   +0.893   0.000  -2.049   0.000];
C2= [ 0.000  -0.351    0.000  +1.646   0.000  -2.587   0.000  -1.678   0.000  -0.170   0.000  -1.003    0.000  -1.306   0.000  -0.492];
% 计算LLR(log-likelihood ratio)
% 认为信道为AWGN信道N(0,1.2)，即SNR = E/N0 = 1/(2*1.2^2)
sigma = 1.2;
Lu  = log((1/sqrt(2*pi*sigma)*exp(-(U -1).^2/(2*sigma*sigma)))./(1/sqrt(2*pi*sigma)*exp(-(U +1).^2/(2*sigma*sigma))));
Lc  = log((1/sqrt(2*pi*sigma)*exp(-(C -1).^2/(2*sigma*sigma)))./(1/sqrt(2*pi*sigma)*exp(-(C +1).^2/(2*sigma*sigma))));
Lc2 = log((1/sqrt(2*pi*sigma)*exp(-(C2-1).^2/(2*sigma*sigma)))./(1/sqrt(2*pi*sigma)*exp(-(C2+1).^2/(2*sigma*sigma))));
 
% 初始化外信息�?0
Lu_ex = zeros(1,length(U));
% 循环迭代
for m =1:16
    Lu1 = Lu + Lu_ex;    
    Lu_ex = Turbo_BCJR(Lu1,Lc,0);
    %Lu_ex  = Lu_new-Lu1;
    Lu2 = Lu + Lu_ex;
    Lu_est = Lu2;
    Lu2 = Turbo_Interleaver(Lu2);
    Lu_ex = Turbo_BCJR(Lu2,Lc2,0);
    %Lu_ex  = Lu_new-Lu2;
    Lu_ex = Turbo_Interleaver(Lu_ex);
    Lu_est = Lu_est + Lu_ex;
    d = sign(Lu_est);
    fprintf('Iteration %d, BER: %d\n',m,sum(abs(d-t))/2);
end
end
 
function [Le] = Turbo_BCJR(Lu, Lc,end_s)
FRAME_LENGTH = length(Lu);
% 计算分支度量（bench metric�?,就是分别计算在时刻i各个输出的概�?
m11 = (Lu+Lc)/2;
m10 = (Lu-Lc)/2;
% 计算前向状�?�度�?(forawrd state metric): 
alpha = zeros(1,(FRAME_LENGTH+2)*4);
alpha(1:4) = -1000*ones(1,4);
alpha(1) = 0;
for k=1:FRAME_LENGTH
    m_t = alpha((k-1)*4+0+1) - m11(k-1+1);
    m_b = alpha((k-1)*4+2+1) + m11(k-1+1);
    alpha(k*4+0+1) = log(exp(m_t)+exp(m_b));    
    m_t = alpha((k-1)*4+0+1) + m11(k-1+1);
    m_b = alpha((k-1)*4+2+1) - m11(k-1+1);
    alpha(k*4+1+1) = log(exp(m_t)+exp(m_b));   
    m_t = alpha((k-1)*4+1+1) - m10(k-1+1);
    m_b = alpha((k-1)*4+3+1) + m10(k-1+1);
    alpha(k*4+2+1) = log(exp(m_t)+exp(m_b));   
    m_t = alpha((k-1)*4+1+1) + m10(k-1+1);
    m_b = alpha((k-1)*4+3+1) - m10(k-1+1);
    alpha(k*4+3+1) = log(exp(m_t)+exp(m_b));   
end
% 计算后向状�?�度量（backward state metris�?
beta = zeros(1,(FRAME_LENGTH+2)*4);
beta((FRAME_LENGTH)*4+0+1:(FRAME_LENGTH)*4+0+4) = -1000*ones(1,4);
beta((FRAME_LENGTH)*4+end_s+1) = 0;
for k=FRAME_LENGTH-1:-1:0
    m_t = beta((k+1)*4+0+1) - m11(k+1);
    m_b = beta((k+1)*4+1+1) + m11(k+1);
    beta(k*4+0+1) = log(exp(m_t)+exp(m_b));    
    m_t = beta((k+1)*4+2+1) - m10(k+1);
    m_b = beta((k+1)*4+3+1) + m10(k+1);
    beta(k*4+1+1) = log(exp(m_t)+exp(m_b));   
    m_t = beta((k+1)*4+0+1) + m11(k+1);
    m_b = beta((k+1)*4+1+1) - m11(k+1);
    beta(k*4+2+1) = log(exp(m_t)+exp(m_b));   
    m_t = beta((k+1)*4+2+1) + m10(k+1);
    m_b = beta((k+1)*4+3+1) - m10(k+1);
    beta(k*4+3+1) = log(exp(m_t)+exp(m_b));   
end
% 计算在输入Y(1:n)的情况下，状态i-1为u',状�?�i为u的联合概�?
Lnew = zeros(1,length(FRAME_LENGTH));
for k = 1:FRAME_LENGTH
    enumerator  = 0;
    denominator = 0;
    t_d = alpha((k-1)*4+0+1) + beta(k*4+0+1) - m11(k-1+1);
    denominator = denominator + exp(t_d);
    t_e = alpha((k-1)*4+0+1) + beta(k*4+1+1) + m11(k-1+1);
    enumerator = enumerator + exp(t_e); 
    t_d = alpha((k-1)*4+1+1) + beta(k*4+2+1) - m10(k-1+1);
    denominator = denominator + exp(t_d);
    t_e = alpha((k-1)*4+1+1) + beta(k*4+3+1) + m10(k-1+1);
    enumerator = enumerator + exp(t_e); 
    t_d = alpha((k-1)*4+2+1) + beta(k*4+1+1) - m11(k-1+1);
    denominator = denominator + exp(t_d);
    t_e = alpha((k-1)*4+2+1) + beta(k*4+0+1) + m11(k-1+1);
    enumerator = enumerator + exp(t_e); 
    t_d = alpha((k-1)*4+3+1) + beta(k*4+3+1) - m10(k-1+1);
    denominator = denominator + exp(t_d);
    t_e = alpha((k-1)*4+3+1) + beta(k*4+2+1) + m10(k-1+1);
    enumerator = enumerator + exp(t_e); 
    Lnew(k) = log(enumerator/denominator);
end
Le = Lnew - Lu;
end


function [out] = Turbo_Interleaver(in)
    out = reshape(reshape(in,4,4)',1,16);
end
