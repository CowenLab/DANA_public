%% RMS
fname = 'amp-D-000.dat';

RHD = INTAN_Read_RHD_file('info.rhd');
sFreq = RHD.frequency_parameters.board_dig_in_sample_rate;

fp = fopen(fname,'rb');
LFP = fread(fp,'int16');
fclose(fp);

figure
plot(LFP)
mid_ix = round(length(LFP)/2);

%%

%  LFP = LFP(mid_ix:end);
% LFP = LFP(1095161:1159677); % forIntan_NewVoltHS_160224_125103\Intan_NewVoltHS_160224_125103
%  LFP = LFP(41716:1459551); % for C:\Cowen\Data\DANA\SSR_160224_160458\SSR_160224_160458
% LFP = LFP(48387:2803225); % C:\Cowen\Data\DANA\Intan_NewVoltHS_160224_125103\Intan_NewVoltHS_160224_125103

LFP = LFP(317972:1728110); % NoSSR_160224_142211\NoSSR_160224_142211

LFP = LFP*RHD.bit_to_uvolt_conversion;
rms(LFP)

