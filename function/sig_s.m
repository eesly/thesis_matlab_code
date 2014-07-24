function [sendS R I]= sig_s(fs)
fca = 294.117e3;        %Âö³åµ÷ÖÆÐÅºÅÔØÆµ
pulse_width = 31*6*2/fca;  
single_width = 6/fca;   

t = 0:1/fs:pulse_width - 1/fs;
s_L = length(t);

t_s = 0:1/fs:single_width - 1/fs;
sig_L = length(t_s);

single_signal = sin(2*pi*fca*t_s);

sendS = zeros(1,s_L);

for num_bit = 1:31
    if(bitget(uint32(hex2dec('b7937045')),num_bit))
      sendS((num_bit-1)*sig_L+1:num_bit*sig_L) = single_signal;
    else
      sendS((num_bit-1)*sig_L+1:num_bit*sig_L) = -single_signal;  
    end   
end
for num_bit = 32:62
    if(bitget(uint32(hex2dec('b7937045')),num_bit-31))
      sendS((num_bit-1)*sig_L+1:num_bit*sig_L) = single_signal;
    else
      sendS((num_bit-1)*sig_L+1:num_bit*sig_L) = -single_signal;  
    end   
end

R = sendS.*cos(2*pi*fca*t)*2;
I = sendS.*sin(-2*pi*fca*t)*2;
% 
% Hd = myfilter(fs,20e3,fs/2);
% R = filter(Hd.Numerator,1,R);
% I = filter(Hd.Numerator,1,I);