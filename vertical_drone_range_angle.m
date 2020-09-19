%%% PROGRAM TO ANALYSE DATA RECORDED USING DCA1000 AND AWR1243 for the file saved in ppro/data1.
clc; clear all; close all;

%% global variables
% Based on sensor configuration.
   numADCBits = 16; % number of ADC bits per sample.
   numADCSamples = 256; % number of ADC samples per chirp.
   numRx = 4; % number of receivers in AWR1243.
   chirpSize = numADCSamples*numRx;
   chirploops= 128; % No. of of chirp loops.  
   numLanes = 2; % do not change. number of lanes is always 4
   isReal = 0; % set to 1 if real only data, 0 if complex data.
   numFrames = 200; 
   numChirps = 1;
   sampleRate = 10; % [Msps]
   timeStep = 1/sampleRate;    % [us]
   chirpPeriod = numADCSamples * timeStep ; % [us]
   plotEnd = numADCSamples * numChirps*numFrames; %for considering all frames.
   Dx = numADCSamples * numChirps ;
   timeEnd = (plotEnd-1) * timeStep;

%% read file
% read .bin file
fid = fopen('adc_data.bin','r');
% adcData = fread(fid, 'int16');
% fclose(fid);
% fileSize = size(adcData, 1);


% % % % % % % % 

adcData = fread(fid, 'int16');
% if 12 or 14 bits ADC per sample compensate for sign extension
if numADCBits ~= 16
    l_max = 2^(numADCBits-1)-1;
    adcData(adcData > l_max) = adcData(adcData > l_max) - 2^numADCBits;
end
fclose(fid);
fileSize = size(adcData, 1);
% real data reshape, filesize = numADCSamples*numChirps
if isReal
    numChirps = fileSize/numADCSamples/numRx;
    LVDS = zeros(1, fileSize);
    %create column for each chirp
    LVDS = reshape(adcData, numADCSamples*numRx, numChirps);
    %each row is data from one chirp
    LVDS = LVDS.';
else
    % for complex data
    % filesize = 2 * numADCSamples*numChirps
    numChirps = fileSize/2/numADCSamples/numRx;
    LVDS = zeros(1, fileSize/2);
    %combine real and imaginary part into complex data
    %read in file: 2I is followed by 2Q
    counter = 1;
    for i=1:4:fileSize-1
        LVDS(1,counter) = adcData(i) + sqrt(-1)*adcData(i+2); 
        LVDS(1,counter+1) = adcData(i+1)+sqrt(-1)*adcData(i+3); 
        counter = counter + 2;
    end
    % create column for each chirp
    LVDS = reshape(LVDS, numADCSamples*numRx, numChirps);
    %each row is data from one chirp
    LVDS = LVDS.';
end

%organize data per RX
adcData = zeros(numRx,numChirps*numADCSamples);
adc_r1=[]; adc_r2=[]; adc_r3=[]; adc_r4=[];

for i=1:256
    adc_r1(:,i) = LVDS(:,i);
    adc_r2(:,i) = LVDS(:,i+256);
    adc_r3(:,i) = LVDS(:,i+512);
    adc_r4(:,i) = LVDS(:,i+768);
end
new_pos = []; new_val = []; new_pos_cal=[];
% easy_positions=zeros(200,200); 
signal_fft_2=[]; complete_signal_fft=[];


% plotting for all frames (Channel 1)
complete_signal_fft_1 = []; complete_signal_fft_2 = [];
complete_signal_fft_3 = []; complete_signal_fft_4 = [];
Doppler_all_1 = []; Doppler_all_2 = []; 
Doppler_all_3 = []; Doppler_all_4 = []; 
% N = 25; 
freq_point = 6; rec1_mean=[];
doppler_axis = linspace(-freq_point,freq_point,chirploops);
range_axis = linspace(-25, 25, 256);
all_fr_1_2=[]; all_fr_2_3=[]; all_fr_1_3=[]; all_fr_2_4=[]; all_fr_1_4=[]; all_fr_3_4=[]; 

for n = 1 : numFrames
    u = n - 1;
    ff_tp = LVDS(1+128*u:128+128*u,:);
    rec1 = ff_tp(:,1:256); rec2 = ff_tp(:,257:512);
    rec3 = ff_tp(:,513:768); rec4 = ff_tp(:,769:1024);
    rec1_fft = fft(rec1'); rec2_fft = fft(rec2'); rec3_fft = fft(rec3'); rec4_fft = fft(rec4');
    fftsrec1 = flip(fftshift(rec1_fft)); fftsrec2 = flip(fftshift(rec2_fft));
    fftsrec3 = flip(fftshift(rec3_fft)); fftsrec4 = flip(fftshift(rec4_fft));
    fftsrec1([1:143],:)=0; fftsrec2([1:143],:)=0; fftsrec3([1:143],:)=0; fftsrec4([1:143],:)=0;
    
    [pos_val_1, pos_1] = max(fftsrec1); [pos_val_2, pos_2] = max(fftsrec2);
    [pos_val_3, pos_3] = max(fftsrec3); [pos_val_4, pos_4] = max(fftsrec4);
    angle_pos_1 = angle(pos_val_1); angle_pos_2 = angle(pos_val_2); 
    angle_pos_3 = angle(pos_val_3); angle_pos_4 = angle(pos_val_4); 
    sin_inv_angle_pos_1_2 = asind(abs(angle_pos_1 - angle_pos_2)/3.14);
    sin_inv_angle_pos_1_3 = asind(abs(angle_pos_1 - angle_pos_3)/3.14);
    sin_inv_angle_pos_1_4 = asind(abs(angle_pos_1 - angle_pos_4)/3.14);
    sin_inv_angle_pos_2_3 = asind(abs(angle_pos_2 - angle_pos_3)/3.14);
    sin_inv_angle_pos_2_4 = asind(abs(angle_pos_2 - angle_pos_4)/3.14);
    sin_inv_angle_pos_3_4 = asind(abs(angle_pos_3 - angle_pos_4)/3.14);
    mean_1_2 = real(mean(sin_inv_angle_pos_1_2)); mean_1_3 = real(mean(sin_inv_angle_pos_1_3));
    mean_1_4 = real(mean(sin_inv_angle_pos_1_4)); mean_2_3 = real(mean(sin_inv_angle_pos_2_3));
    mean_2_4 = real(mean(sin_inv_angle_pos_2_4)); mean_3_4 = real(mean(sin_inv_angle_pos_3_4));
    
    all_fr_1_2 = [all_fr_1_2 mean_1_2]; all_fr_1_3 = [all_fr_1_3 mean_1_3];
    all_fr_1_4 = [all_fr_1_4 mean_1_4]; all_fr_2_3 = [all_fr_2_3 mean_2_3];
    all_fr_2_4 = [all_fr_2_4 mean_2_4]; all_fr_3_4 = [all_fr_3_4 mean_3_4];
      
    angle1 = angle(fftsrec1); angle2 = angle(fftsrec2); angle3 = angle(fftsrec3); angle4 = angle(fftsrec4); 
    angle_diff_1_2 = abs(angle1-angle2); sin_inv_angle_1_2=asind(angle_diff_1_2/3.14);
    angle_diff_1_3 = abs(angle1-angle3); sin_inv_angle_1_3=asind(angle_diff_1_3/(2*3.14));
    angle_diff_1_4 = abs(angle1-angle4); sin_inv_angle_1_4=asind(angle_diff_1_4/(3*3.14));
    angle_diff_2_3 = abs(angle2-angle3); sin_inv_angle_2_3=asind(angle_diff_2_3/3.14);
    angle_diff_2_4 = abs(angle2-angle4); sin_inv_angle_2_4=asind(angle_diff_2_4/(2*3.14));
    angle_diff_3_4 = abs(angle3-angle4); sin_inv_angle_3_4=asind(angle_diff_3_4/3.14);

    figure(1); plot(range_axis,abs(fftsrec1(:,128))); grid on;
    figure(2); plot(range_axis,abs(fftsrec2(:,128))); grid on;
%     drawnow;
end
final_angle_1_2 = mean(all_fr_1_2); final_angle_1_3 = mean(all_fr_1_3);
final_angle_1_4 = mean(all_fr_1_4); final_angle_2_3 = mean(all_fr_2_3);
final_angle_2_4 = mean(all_fr_2_4); final_angle_3_4 = mean(all_fr_3_4);
final_angle_avg = (final_angle_1_2+final_angle_1_4+final_angle_2_4+final_angle_1_3+final_angle_2_3+final_angle_2_3)/6;
final_angle_avg_on_3 = (final_angle_1_2+final_angle_2_3+final_angle_3_4)/3

% 
% 
% % ========= (3) MUSIC ALGORITHM ========= %
%     % Sample covariance matrix
%     Rxx = x*x'/L;
%     d = 0.5;% Distance between elements in Wavelengyh
%     [Q ,D] = eig(Rxx); %Compute eigendecomposition of covariance matrix
%     [D, I] = sort(diag(D),1,'descend'); %Find r largest eigenvalues
%     Q = Q(:,I); %Sort the eigenvectors to put signal eigenvectors first
%     Qs = Q(:,1:M); %Get the signal eigenvectors
%     Qn = Q(:,M+1:Ant); %Get the noise eigenvectors
%     % MUSIC algorithm
%     % Define angles at which MUSIC “spectrum” will be computed
%     angles = (-28:1:28);
%     for k = 1 : length(angles)
%         a1(:,k) = exp(-1i*2*pi*(d*(0:Ant-1)'*sin([angles(k).']*pi/180)));
%     end
%     for k = 1 : length(angles)
%         %Compute MUSIC “spectrum”
%         music_spectrum(k) = (a1(:,k)'*a1(:,k))/(a1(:,k)'*(Qn*Qn')*a1(:,k));
%     end
%     plot(angles,abs(music_spectrum));
%     grid on; title('MUSIC Spectrum Angle estimation'); xlabel('Angle in degrees'); 
%     hold on;
% 
% 
% 
% 
% 




% s_fft = mean(signal_fft_2');
% new_pos_cal(isnan(new_pos_cal))=0;
% % new_pos_cal_2 = mean(new_pos_cal);
% new_pos_cal_2 = max(new_pos_cal);
% [mxpp_new, lopp_new] = max(new_pos_cal_2);
% % subplot(2,2,4);
% plot(range,new_pos_cal_2); 
% plot(range(lopp_new),mxpp_new,'r--o','MarkerSize',10);
% title(sprintf('Average of all frame, Drone is at %f', range(lopp_new)));
% axis tight; grid on; xlabel('Range (m)'); ylabel('Amplitude');
% grid minor; hold on;

% 
% %%% Micro Doppler Signatures
% frameTime = 0.040;   % 40ms per frame
% freq_point = 6; 
% totalTime = frameTime*numFrames;
% 
% n1 = 80; n2 = n1+1; n3 = n1+2; n4 = n1+3;
% md1 = complete_signal_fft(n1,:); md2 = complete_signal_fft(n2,:); 
% md3 = complete_signal_fft(n3,:); md4 = complete_signal_fft(n4,:); 
% new_md1 = []; new_md2 = []; new_md3 = []; new_md4 = []; 
% delta = 12; numb = 100; i_end = (length(md1) - numb)/delta + 1;
% for i = 1 : i_end
%     new_md1(:,i) = fft(md1( 1 + (i-1)*delta : (i-1)*delta + numb ));
%     new_md2(:,i) = fft(md2( 1 + (i-1)*delta : (i-1)*delta + numb ));
%     new_md3(:,i) = fft(md3( 1 + (i-1)*delta : (i-1)*delta + numb ));
%     new_md4(:,i) = fft(md4( 1 + (i-1)*delta : (i-1)*delta + numb ));    
% end
% 
% figure('Name','micro doppler');
% time_axis = linspace(0,totalTime,i_end); freq_axis = linspace(-freq_point,freq_point,numb);
% subplot(2,2,2); imagesc(time_axis,freq_axis,mag2db(abs(fftshift(new_md1))));
% xlabel('Time (s)'), ylabel('Frequency (Hz)'); axis tight; grid on; grid minor;
% % title(['  drone at ' num2str(distance(n1))],'Interpreter','tex');