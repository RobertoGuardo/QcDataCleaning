% SCRIPT AIM: Datasete analysis and cleaning before using into the Qc code 
%
% Author:  Roberto Guardo
%
%           Apr 2018
%           
%           Last Update: February 2020
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Folders connection and import
disp('START')
clear
clc 
% warning off all
addpath('Utilities_Matlab/');
addpath('Traces');

load tP_20283.txt; % Time picking P
tP = tP_20283;

% End begin


%% STEP 0 - Setup

% AREA: Decepcion Island ( DCP )

nW = 1.5;   % noise Windows
cf = 3;     % central frequency (from Prudencio et al., 2015)
cW = 4;     % coda Windows
mL = 10;    % total length of the signal (maxLength) 
sn = 0;     % start time to consider noise - before P-arrival
 
srate = 100; % sampling rate (time(n)-time(n-1))^-1;

% Iteration, choose from 0 to 4
% 0 = Full dataset 20283 events
% 1 = 14972 rays
% 2 = 13105 rays
% 3 = 7895 rays
% 4 = 7196 rays
it = 0;
% End of Setup0

traceIt0 = 'traces_it0_20283.txt';


%% Check single trace

switch it
    case 1
        lista = textread('traces_it1_14972.txt','%q','delimiter','\n');
    case 2
        lista = textread('traces_it2_13105.txt','%q','delimiter','\n');
    case 3
        lista = textread('traces_it3_7895.txt','%q','delimiter','\n');
    case 4
        lista = textread('traces_it4_7196.txt','%q','delimiter','\n');
    otherwise
        %lista = textread('traces_it0_20283.txt','%q','delimiter','\n');
        lista = textread(traceIt0,'%q','delimiter','\n');
end

ll = length(lista);

i = 1; %1870; %1931; %1933; %1939; %1944;
listex = lista{i,1};    % reading seismograms
    
[trace] = textread(listex,'%f');  % get 2-columns files
time = trace(1:2:end-1);    % first column is time
ampl = trace(2:2:end);      % second is amplitude
    
time = time(1:srate*mL);    % cut the trace according to the tW
ampl = ampl(1:srate*mL);

posamp = abs(ampl);
spAmp = fft(ampl);
spPosamp = abs(spAmp);

% ------- Visualize Trace, Trace with abs values and Fft of the trace
    figure('Name',char(lista(i)));
    h1 = subplot(3,1,1); 
    subplot(3,1,1);
    plot(time,ampl);
    title(['Trace number: ', num2str(i)]);
    axes(h1);
    line([tP(i,2) tP(i,2)],get(h1,'YLim'),'Color',[1 0 0]);
    subplot(3,1,2);
    plot(time,posamp);
    title('Positive Amplitude');
    subplot(3,1,3);
    plot(time,spPosamp);
    title('FFT Positive Amplitude');
% ------------------------------------------------------------------- 


% ------- Visualize Trace, its spectrogram and Fft of the trace
    figure('Name', ['Spectrogram ',char(lista(i)), ' ID: ', char(num2str(i))]);
    subplot(3,1,1);
    plot(time,ampl);
    title(['Trace number: ', num2str(i)]);
    subplot(3,1,2);
    spectrogram(ampl,[],[],[],srate,'yaxis');
    title(['Spectrogram: ', num2str(i)]);
    subplot(3,1,3);
    plot(time,spPosamp);
    title('FFT Positive Amplitude');
% -------------------------------------------------------------------


% ------- Visualize Trace, its spectrogram and Fft of the trace
    f = (0:length(ampl)-1)*(srate/length(ampl));
    power = abs(spAmp).^2/length(ampl);
    subplot(3,1,1);
    spectrogram(ampl,[],[],[],srate,'yaxis');
    title(['Spectrogram: ', num2str(i)]);
    subplot(3,1,2);
    plot(f(1:numel(f)/2),power(1:numel(power)/2));
    xlabel('Frequency');
    ylabel('Power');
    title(['PowerSpectrum: ', num2str(i)]);
    subplot(3,1,3);
    loglog(f(1:numel(f)/2),power(1:numel(power)/2));
    xlabel('Frequency');
    ylabel('Power');
    title('Log/Log PowerSpectrum');
    
    [maxPow,iMaxPow] = max(power);
    maxPowF = f(iMaxPow);
% -------------------------------------------------------------------

% end CHECK SINGLE TRACE


%% Visualise traces 
% (Move it after "Trace Analysis")

switch it
    case 1
        lista = textread('traces_it1_14972.txt','%q','delimiter','\n');
    case 2
        lista = textread('traces_it2_13105.txt','%q','delimiter','\n');
    case 3
        lista = textread('traces_it3_7895.txt','%q','delimiter','\n');
    case 4
        lista = textread('traces_it4_7196.txt','%q','delimiter','\n');
    otherwise
        lista = textread('traces_it0_20283.txt','%q','delimiter','\n');
end
ll = length(lista);
% load('checkTrace.mat');
i=1;
j=1;
goodTraces = 0;
for i = 1:ll
    
    listex = lista{i,1};    % reading seismograms
    
    [trace] = textread(listex,'%f');  % get 2-columns files
    time = trace(1:2:end-1);    % first column is time
    ampl = trace(2:2:end);      % second is amplitude
    
    time = time(1:srate*mL);    % cut the trace according to the tW
    ampl = ampl(1:srate*mL);
    

%    charS = checkTraceStr(i,16);
%    charSt = char(charS);

%     % IT 1 --- Zeros    
%     if checkTrace(i,2) > 40                         
%         continue
%     % IT 2 --- Coda Window shorter than 4 sec
%     elseif infoPosamp(i,3) == 0                     
%         continue
%     % IT 3 --- Coda window Correlation Coefficients
%     elseif checkTrace(i,12) <= 0.4                  
%         continue
%     % IT 4 --- Mean Aplitude (spikes)
%     elseif checkTrace(i,10) >= 94                    
%         continue
%     % IT 5 --- Bad Stations
%     elseif (charSt(1) == "W") || (charSt(1:3) == "E09") 
%          continue
%     end
    
%   goodTraces = goodTraces + 1;
  
%     h(j) = subplot(6,3,j);
%     subplot(6,3,j);
%     f = plot(time,ampl);
%     title([lista(i),i]);
%     axes(h(j));
%     line([tempi(i,2) tempi(i,2)],get(h(j),'YLim'),'Color',[1 0 0]);
  
    subplot(6,3,j);
    f = plot(time,ampl);
    title(['Traza número: ', num2str(i), '- Estación: ',listex(12:15)]);
    xlabel('Tiempo (s)');
    ylabel('Amplitud');
   %title(listex(12:15));
   %title([lista(i),i])
   %title(i);

    %if listex(12:15) == 'K11N'                
        j = j + 1;
        if j <= 18
            continue
        else
            w = waitforbuttonpress;
            if w == 0
                disp('Button click')
            else
                disp('Key press')
            end
            j = 1;
        end
    %end
end

% end VISUALISE TRACES


%% Trace analysis

%%%%%%%%% Here we find the parameters to discriminate 
%%%%%%%%% whether a trace is good or not

display('Start - Trace analysis');

lista = textread('traces_it0_20283.txt','%q','delimiter','\n');
ll = length(lista);    
    
corrCoef = zeros(ll,6);
checkTrace = zeros(ll,13);
staLta = zeros(ll,7);

cf = 3;             % central frequency, the first value is equal to cf+3
threshold = 0.4;    % Minimum RZZ
index = 0;
keepCf = cf;
coda = 4;           % Select 1,2,3 or 4

for i = 1:ll
    index = index+1;
    display(index);
    
    listex = lista{i,1};              % reading seismograms
    [trace] = textread(listex,'%f');  % get 2-columns files
    
    time = trace(1:2:end-1);    % first column is time
    ampl = trace(2:2:end);      % second is amplitude
    
    time = time(1:srate*mL);    % trace cutting according to the tW
    ampl = ampl(1:srate*mL);    % trace cutting according to the tW
    
    posamp = abs(ampl);
    
    f = (0:length(ampl)-1)*(srate/length(ampl));
    power = abs(fft(ampl)).^2/length(ampl);
    [maxPow,iMaxPow] = max(power);
    maxPowF = f(iMaxPow);
    
    cf = keepCf;        %Necessary to reset the cf before the for loop
    nf=[4 4 4 8 8 8];

    for n = 1:6
        
        cf = cf + 3;                        % Increase frequency
        Wn = ([cf-cf/3 cf+cf/3]/50);        % Frequency band
        
        [z,p,k] = butter(6,Wn,'bandpass');  % Butter filter
        [sos,g] = zp2sos(z,p,k);            % Convert to SOS form
        Hd = dfilt.df2tsos(sos,g);          % Create a dfilt object   
    
               
        % NOISE %---------------------------------------------------------%
        intn = nW*srate;                    % Number of sample for noise
		cursorN1 = 1;                       % Starting sample for noise
        cursorN2 = cursorN1 + intn-1; 		% End-sample of noise-window
        tsisman = ampl(cursorN1:cursorN2);	% Tapering noise
        
        % ENVELOPE %------------------------------------------------------%
        fsisman = filter(Hd,tsisman);		% filtering noise
        hspn = hilbert(fsisman);			% hilbert transform
        spmn = abs(hspn.^2);				% hilbert
        spsn = smooth(spmn,nf(n)/cf*srate); % smooth hilbert  
        
        %ENTIRE - We now compute the entire envelope energy %-------------%
		fsismac = filter(Hd,ampl);          % filter coda
        hsp = hilbert(fsismac);             % hilbert
        spm = abs(hsp.^2);                  % ms of the filtered waveform
        sps = smooth(spm,nf(n)/cf*srate);   % smooth ms 
        
        % CODA %----------------------------------------------------------% 
        
        switch coda
            case 1 % Coda from 4s to the end of the trace
                sps2 = sps(4*srate:length(ampl));
                
            case 2 % Coda from 6s to the end of the trace
                sps2 = sps(6*srate:length(ampl)); 

            case 3 % Coda from 4s to 8s
                sps2 = sps(4*srate:4*srate*2);
 
            otherwise % Coda from tS to (tS + 4s)
                [maxsps,imaxsps] = max(sps);
                 cLim = imaxsps + cW*srate;          % coda limit
                
                % NB: cLim is previous the "if" in order to get its last
                % value, even if it is > length(ampl). 
                % The next "if" will give 0 to the relative corrcoeff.

                if (cLim - length(ampl)) >= 40
                    continue
                else
                    cLim = length(ampl); % Cut the coda if exceeds
                end

                sps2 = sps(imaxsps:cLim);
                
        end
        
		% CORRELATION COEFFICIENTS %--------------------------------------% 
        lspm = length(sps2);
        tm = (imaxsps:lspm+imaxsps-1)'/srate;
        lspm3 = log(sps2.*tm.^1.5)/2/pi/cf; % linearisation
                
        Rz=corrcoef(tm,lspm3);
        if Rz == 1
            Rz(1,2) = 1;
        end
        
        corrCoef(i,n) = abs(Rz(1,2)); 
        
    end % ----------------------------------------------------------------------------------------------------------
    
    if cLim > length(ampl)
        cLim = length(ampl);
    end
    

    % Define sample trace     
    if i == 1 % (check before choose the sample)
        sample = posamp;
    end
    
    verXcorr = xcorr(sample,posamp);
    
%     % IT 2 
%     [maxPosamp,imaxPosamp] = max(posamp);
%     if imaxPosamp > 600
%         continue
%     end
    
    checkTrace(i,1) = i; % ID
    checkTrace(i,2) = sum(posamp<=0.5)/length(posamp)*srate; % zeros
    checkTrace(i,3) = mean(posamp(1:nW*srate-1)); % nW amplitude mean
    checkTrace(i,4) = mean(posamp((nW*srate+1):((nW*srate+1)+3*srate))); % signal
    checkTrace(i,5) = checkTrace(i,4)/checkTrace(i,3); % Signal to Noise ratio
    checkTrace(i,6) = mean(posamp(imaxsps:cLim)); % mean of the Coda
    checkTrace(i,7) = maxPow; % Max PowerAmplitude Mean Jumps values
    checkTrace(i,8) = maxPowF; % Frequency at MaxPow
    checkTrace(i,9) = mean(verXcorr); % mean of the xCorr values
    checkTrace(i,10) = sum(posamp>500); % sum of the amplitude value >500 (spikes)
    for j = 1:length(posamp)-1
        jumpMean(j) = mean(abs(posamp(j+1)-posamp(j))); 
        checkTrace(i,11) = mean(jumpMean); % mean of the jumps values
    end
    checkTrace(i,12) = mean(corrCoef(i,:));
    
end

% STA/LTA (Short Time Average over Long Time Average) %-------------------%

ltW = (numel(posamp)+1-1*srate);    % Long Time Window

for i = 1:ll
   
    listex = lista{i,1};                % reading seismograms
    [trace] = textread(listex,'%f');    % get 2-columns files
    
    ampl = trace(2:2:end);      
    ampl = ampl(1:srate*mL);
    
    posamp = abs(ampl);
    
    [maxPosamp,imaxPosamp] = max(posamp);
    infoPosamp(i,1) = i;
    infoPosamp(i,2) = imaxPosamp;
    infoPosamp(i,3) = (infoPosamp(i,2) <= 600); % maxPosAmp > 600 -> 0 
                                                % coda less than 4 s

stW = 1;

    for s = 1:6

        staLta(i,s) = mean(posamp(stW:(stW+nW*srate)))/mean(posamp(1:ltW));
        stW = stW + nW*srate; % Moving the short time windows with a step equal to nW

    end
    
    staLta(i,7) = std(staLta(i,[1:6]));
    checkTrace(i,13) = staLta(i,7);
    
end

for k = 1:numel(checkTrace(1,:))-1
    abouTraces(k) = mean(checkTrace(:,k+1));
end

% Add the station name %--------------------------------------------------%

checkTraceStr = string(checkTrace); % Create a checkTrace string array

for i = 1:ll
    listex = lista{i,1};
    checkTraceStr(i,14) = listex(12:15);
end


% for i = 1:ll
%     if str2num(char(checkTraceStr(i,16)))<threshold
%         blwThresh(i,1) = checkTraceStr(i,16);
%         blwThresh(i,2) = checkTraceStr(i,17);
%     else
%         %continue
%         blwThresh(i,2) = "ok";
%     end    
% end
% 
% % find unique staz value and count them
% % then compare with the staz list 
% 
% blwThresh = rmmissing(blwThresh);
% [stazNameBT,ia,ic] = unique(blwThresh(:,2));
% name_counts = accumarray(ic,1);
% BT_stat_counts = [stazNameBT, name_counts]; %which stations and how many with value below trheshold

% check the time of the shot and if the trace's quality
% find a relation to discard bad traces

save checkTrace.mat


% end of Trace evaluation

    
%% New dataset creation

lista = textread('traces_20283.txt','%q','delimiter','\n');
ll = length(lista);
listaR = textread('rays_20283.txt','%q','delimiter','\n');
llR = length(listaR);

tot = 1;

i=1;
for i = 1:ll
    
    listex = lista{i,1};    % reading seismograms
    listexR = listaR{i,1};  % reeading traces
    
    
    % Memorandum 
    % checkTrace(1) = ID
    % checkTrace(2) = Zeros
    % checkTrace(3) = Mean Noise Window 
    % checkTrace(4) = Mean Signal 
    % checkTrace(5) = Mean Signal to Noise 
    % checkTrace(6) = Mean Coda windows
    % checkTrace(7) = Max PowerAmplitude
    % checkTrace(8) = Frequency at MaxPow
    % checkTrace(9) = xCorr mean values
    % checkTrace(10) = Sum of the amplitude value > 500 (spikes)
    % checkTrace(11) = Mean of the jumps values
    % checkTrace(12) = Mean CorrelationCoefficents
    % checkTrace(13) = Standard Devation Sta/Lta windows

     
    charS = checkTraceStr(i,14);
    charSt = char(charS);

        
    if checkTrace(i,2) > 40                         % IT 1 --- Zeros
        continue
    elseif infoPosamp(i,3) == 0                     % IT 2 --- CodaW < 4 s
        continue
    elseif checkTrace(i,12) <= 0.6                  % IT 3 --- CorrCoeff
        continue
%     elseif checkTrace(i,10) >= 94                   % IT x --- Mean Spikes
%         continue
    elseif (charSt(1) == "W") || (charSt(1) == "F") % IT 4 --- Bad Stations
         continue
    end

     
    fileOut = fopen('traces_it4.txt','a');
    fprintf(fileOut,'%s\n',listex);
    fclose(fileOut);
    
    fileOut = fopen('rays_it4.txt','a');
    fprintf(fileOut,'%s\n',listexR);
    fclose(fileOut);
    
    tot = tot+1;
end

% List "Station-Events"
% [statName,SNia,SNic] = unique(checkTraceStr(:,15));
% statNameCount = accumarray(SNic,1);
% infoStat = [statName, statNameCount];



% end Creation of clean dataset


