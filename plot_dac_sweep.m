%Modified version of Cheng's program radioactive_source_spectrum_20111202.m

clear
overflow_compensate=true; %Assume that each time overflow has occured, it has happened only once and only in the lowest bin. Use with care!
overflow_bin=1; % The bin we add counts to if overflow occurs. Only used if overflow_compensate is true.
store_all_counter_values=true; %Set to false to only store the sum of the counters in each channel instead of the individual frames (saves memory).
frames = 14;
energybins = 8; % represent 8 counters.
dac_values = 0:1:20; %Remember to change this depending on which sweep you want to load
scanned_dac_nr = 0; %Note: starts counting on 0
skipped_frames_in_beginning=4;
repetitions = 1;
ch =1:160; %Note that ch starts numbering at 1
channels=length(ch);
unconnected_channels = [];
connected_channels = setdiff(ch,unconnected_channels);

%We may want to keep track of the folder one level above the raw data on if we wnat to save plots and don't want them mixed with the raw data files.
%scan_folder = './scan/'; %Should end with a slash
%rootpath1 = [scan_folder 'raw_data/'];

rootpath1 = './test_data/';

total_counts = zeros(length(dac_values), channels,energybins); %Stores the sum of the counts in different frames and repetitions
if(store_all_counter_values==true)
    all_counter_values=zeros(length(dac_values),channels,frames,repetitions,energybins);
end

%load via_channels via_channel_numbers %Only if we want to separate via channels and nonvia channels when plotting

asics_to_plot=[0 1 2 3 4]; %Note: starts counting at 0
%Note that the count matrices are overwritten when the program goes on to
%the next asic.

for entry_in_plot_list_nr = 1:length(asics_to_plot)
    asic_nr=asics_to_plot(entry_in_plot_list_nr);

    %via_channels_current_asic = via_channel_numbers((via_channel_numbers>160*asic_nr) & (via_channel_numbers<=160*(asic_nr+1)))-160*asic_nr;
    %Note: here channel numbering starts with 1.
    %nonvia_channels_current_asic =setdiff(ch,via_channels_current_asic);
    
    for entry_in_dac_list_nr = 1:length(dac_values)
        DAC = dac_values(entry_in_dac_list_nr);
        disp(DAC);

        for repetition_nr = 1:repetitions

            if DAC < 10
                path2 = 'DAC0';
            else
                path2 = 'DAC';
            end

            path3 = num2str(DAC); %%%% change to correonding kVp spectrum for imaging
            if repetition_nr <= 10
                path4 = 'Repeat0';
            else
                path4 = 'Repeat';
            end
            path5 = num2str(repetition_nr-1);

            name = [rootpath1,path2,path3,path4,path5,sprintf('DAC%dScanSelected_asic%d.asic',scanned_dac_nr,asic_nr)];

            [counters, measurements_with_overflow, ~] = loadFramesFromFile(name, channels, frames, energybins);
            %counters has size: (channels, frames,energybins)
            %measurements_with_overflow has size: (channels, frames)
            
            if(overflow_compensate==true)
                counters_in_overflow_bin = counters(:,:,overflow_bin);
                counters_in_overflow_bin(measurements_with_overflow==1)=counters_in_overflow_bin(measurements_with_overflow==1)+255; %Note that the LFSR counters only have 255 states, not 256.
                counters(:,:,overflow_bin) = counters_in_overflow_bin;
            end
            if(store_all_counter_values==true)
                all_counter_values(entry_in_dac_list_nr,:,:,repetition_nr,:) = counters;
            end
            total_counts(entry_in_dac_list_nr,:,:) = permute(total_counts(entry_in_dac_list_nr,:,:),[2 3 1]) + permute(sum(counters,2),[1 3 2]);
        end
    end

    figure;
    axes('fontsize',12);
    plot(dac_values,sum(sum(all_counter_values(:,:,skipped_frames_in_beginning+1:end,:,scanned_dac_nr+1),3),4)); %To plot a single bin
    %plot(dac_values,sum(sum(sum(all_counter_values(:,:,skipped_frames_in_beginning+1:end,:,:),3),4),5)); %To plot the sum over all bins

    %plot(dac_values,sum(sum(all_counter_values(:,via_channels_current_asic,skipped_frames_in_beginning+1:end,:,scanned_dac_nr+1),3),4),'r'); hold all %To plot a single bin, only via channels
    %plot(dac_values,sum(sum(all_counter_values(:,nonvia_channels_current_asic,skipped_frames_in_beginning+1:end,:,scanned_dac_nr+1),3),4),'b'); %To plot a single bin, only nonvia channels
    xlabel('DAC value');
    ylabel('Counts');
    title(sprintf('ASIC%d',asic_nr));
    grid on
    %saveas(gcf,sprintf([scan_folder 'asic%d.fig'],asic_nr))
end
