frq_characteracter = getDefaultfrq_characteracter();
% Initialize GUI
main_gui = createGUI(frq_characteracter);
% Create the GUI 
function main_gui = createGUI(frq_characteracter)
    % Create a figure.
    main_gui = figure('Position', [100, 100, 500, 300], 'Units', 'pixels');
    % field for entering the in_sentence.
    InputField = uicontrol('Style', 'edit', 'String', '', 'Units', 'pixels', 'Position', [20, 200, 400, 20]);
    % button for encoding the in_sentence.
    Button_Encode = uicontrol('Style', 'pushbutton', 'String', 'Encode and Play', 'Units', 'pixels', 'Position', [20, 170, 150, 20]);
    % button saving the file(.wav )
    Button_save = uicontrol('Style', 'pushbutton', 'String', 'Save as .wav', 'Units', 'pixels', 'Position', [200, 170, 150, 20]);
    % button for uploading an audio file.
    Button_upload = uicontrol('Style', 'pushbutton', 'String', 'Upload Audio File', 'Units', 'pixels', 'Position', [20, 140, 150, 20]);
    % button for decoding using fft.
    Button_decode_fft = uicontrol('Style', 'pushbutton', 'String', 'Decode (Frequency Analysis)', 'Units', 'pixels', 'Position', [200, 140, 200, 20]);
    % button for decoding with bandpass filtering.
    Button_decode_bandpass = uicontrol('Style', 'pushbutton', 'String', 'Decode (Bandpass Filtering)', 'Units', 'pixels', 'Position', [200, 110, 200, 20]);
    % Add all button on the gui  
    set(Button_Encode, 'Callback', {@encode_inputt, main_gui, frq_characteracter});
    set(Button_save, 'Callback', {@save_file_wav, main_gui, frq_characteracter});
    set(Button_upload, 'Callback', {@upload_audiofile, main_gui});
    set(Button_decode_fft, 'Callback', {@decode_fft_button, main_gui});
    set(Button_decode_bandpass, 'Callback', {@decode_filter_bandpass, main_gui});
set(Button_decode_fft, 'Callback', {@decode_fft_button, main_gui, frq_characteracter});
    guidata(main_gui, struct('InputField', InputField, 'Button_Encode', Button_Encode, 'Button_save', Button_save, ...
        'Button_upload', Button_upload, 'Button_decode_fft', Button_decode_fft, 'Button_decode_bandpass', Button_decode_bandpass));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%function of button 
% button for decoding with fft
function decode_fft_button(~, ~, main_gui, frq_characteracter)
    % Get sentence from (path file audio) field and buttons from the figure window.
    data = guidata(main_gui);
    InputField = data.InputField;
    % take path are selected
    audio_file = get(InputField, 'String');
    % Decode using fft
    decodedinput = decode_in_signal(audio_file, frq_characteracter);
    % show output 
    msgbox(['Decoded String (Frequency Analysis): ' decodedinput], 'Decoding Result');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% button saves the audio file as a .wav file.
function save_file_wav(~, ~, main_gui, frq_characteracter)
    % Get the edit field and buttons from the figure window.
    data = guidata(main_gui);
    InputField = data.InputField;
    % Get the in_sentence 
    in_sentence = get(InputField, 'String');
    % Encode the in_sentence.
    encodedin_signal = encode_in_sentence(in_sentence, frq_characteracter);
    try
        % Ask the user to select a name of file &  path
        [filename, pathname] = uiputfile('*.wav', 'Save as .wav');
        % Check if the user canceled the operation or no
        if isequal(filename, 0) || isequal(pathname, 0)
            disp('Save operation canceled by user.');
            return;
        end
        % Concatenate the full path
        path = fullfile(pathname, filename);
        % Save the encoded in_signal as a .wav file.
        fs = 8000;
        audiowrite(path, cell2mat(encodedin_signal), fs);
        % Convert text to speech and play
        player = audioplayer(cell2mat(encodedin_signal), fs);
        play(player);
        [y, Fs] = audioread(path);
        sound(y, Fs);
    catch
        disp('Error saving the .wav file.');
    end
end
% button for encode
function encode_inputt(~, ~, main_gui, frq_characteracter)
    % Get sentence from  field and buttons from the figure window.
    data = guidata(main_gui);
    InputField = data.InputField;
    % Get the in_sentence from field.
    in_sentence = get(InputField, 'String');
    % Encode input.
    encodedin_signal = encode_in_sentence(in_sentence, frq_characteracter);
    % Convert the in_signal to an audio& play .
    in_signal_attime = cell2mat(encodedin_signal);
    sound(in_signal_attime, 8000);
player = audioplayer(in_signal_attime, 8000);
play(player);
end
% button for uploading an audio file.
function upload_audiofile(~, ~, main_gui)
    % Get input from field and buttons
    data = guidata(main_gui);
    InputField = data.InputField;
    % choose an audio file want to upload.
    [filename, pathname] = uigetfile('*.wav', 'Select an audio file');
    % Check if  canceled the operation if yes show masseg canceled 
    if isequal(filename, 0) || isequal(pathname, 0)
        disp('File selection canceled by user.');
        return;
    end
 % show the selected file path field.
    set(InputField, 'String', fullfile(pathname, filename));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% button for decoding with bandpass 
function decode_filter_bandpass(~, ~, main_gui)
    % Get sentance from  field and buttons from the main window.
    data = guidata(main_gui);
    InputField = data.InputField;
    % Get the path of selected audio file.
    audio_file = get(InputField, 'String');
    % Decode using bandpass 
    decodedinput = decode_filtering(audio_file);
    % Create a new GUI window & show output.
    guiFigure = figure('Name', 'Decoded String GUI', 'NumberTitle', 'off', 'Position', [100, 100, 300, 100]);
    uicontrol('Style', 'text', 'Position', [10, 40, 280, 30], 'String', ['Decoded String: ' decodedinput]);
    % Save the GUI figure handle in the 'UserData' 
    data.guiFigure = guiFigure;
    guidata(main_gui, data);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%function to encode & decode
% Function to encode a in_sentence
function encodedin_signal = encode_in_sentence(in_sentence, frq_characteracter)
    encodedin_signal = {};
    % Convert the in_sentence to lowercase
    in_sentence = lower(in_sentence);
    for character = in_sentence
        % Check if the character consist in frq_characteracter
        if isfield(frq_characteracter, character)
            % define numofsample ,fs,frequcy for each char
            frq = frq_characteracter.(character);
           fs=8000;
           m=0:320;        
% Generate the in_signal using the cos
in_signal = cos((2*pi*frq(1)*m)/fs) + cos((2*pi*frq(2)*m)/fs) + cos((2*pi*frq(3)*m)/fs);
            encodedin_signal{end+1} = in_signal;        
       elseif character == ' ' % Check if char== space         
            frq = frq_characteracter.space;      
% Generate the in_signal 
in_signal = cos((2*pi*frq(1)*m)/fs) + cos((2*pi*frq(2)*m)/fs) + cos((2*pi*frq(3)*m)/fs);
            % add the encoded in_signal to the output
            encodedin_signal{end+1} = in_signal;
        else
           % dont care for anther charcter not in frq_characteracter 
    errorMsg = ['characteracter ''' character ''' not supported.'];
    errordlg(errorMsg, 'Error', 'modal');
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  decode using narrow bandpass filters 
function decodedinput = decode_filtering(audio_file, frq_characteracter)
    % Check if frq_characteracter is not provided
    if nargin < 2
        frq_characteracter = getDefaultfrq_characteracter();
    end
    % Read the audio file
    [audio, fs] = audioread(audio_file);
    decodedinput = '';
  % Frame size
    windowize = 320;  
    % div audio to window (size=320)
    for i = 1:windowize:length(audio)-windowize
        audioFrame = audio(i:i+windowize-1);
        outputs = zeros(length(fieldnames(frq_characteracter)), 1);
        % Apply narrow bandpass 
        y = 1;
        for characteracter = fieldnames(frq_characteracter)'
            characterFrequencies = frq_characteracter.(characteracter{1});
            filters = cell(1, numel(characterFrequencies));
% Create narrow bandpass filter (iir)for each frequency (low,mid,high)
for k = 1:numel(characterFrequencies)
    filters{k} = designfilt('bandpassiir', 'FilterOrder', 10, ...
        'HalfPowerFrequency1', max(characterFrequencies(k) - 10, 1), 'HalfPowerFrequency2', min(characterFrequencies(k) + 10, fs/2), ...
        'SampleRate', fs);
end
            % Apply each filter to the current frame & calculate output
            filteredOutputs = zeros(1, numel(filters));
            for k = 1:numel(filters)
                filteredFrame = filter(filters{k}, audioFrame);
                filteredOutputs(k) = sum(abs(filteredFrame));
            end
            outputs(y) = sum(filteredOutputs);
            y = y + 1;
        end
        % Find max output & convert index of character to frequncy
        [~, characterIndex] = max(outputs);
        characterr = char(96 + characterIndex); 
        if (characterr == '{')
        characterr = 'y';
        end
        decodedinput = [decodedinput, characterr];
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function decoded_string = decode_in_signal(audio_file, frq_characteracter)
    % Load the audio in_signal
    [audio_data, Fs] = audioread(audio_file);
    display(Fs);
    % Define frame parameters
    duration = 0.04; % 40 milliseconds
    window_size = round(duration * Fs);
    display(window_size);
    number_window = floor(length(audio_data) / window_size);
    % Initialize variables
    window = zeros(window_size, number_window);
    decoded_characteracters = cell(1, number_window);
    % Reshape the audio data into window
    for frame_index = 1:number_window
        window(:, frame_index) = audio_data((frame_index - 1) * window_size + 1 : frame_index * window_size);
    end
    % Decode each frame
    for frame_index = 1:number_window
        frame = window(:, frame_index);
        % Apply Fourier transform
        fft_output = fft(frame(1:320));
        fft_mag = abs(fft_output);
        % Find the three highest peaks (excluding DC component)
        [~, peak_indices] = findpeaks(fft_mag(2:150), 'SortStr', 'descend', 'NPeaks', 3);
        % Convert peak indices to frequencies
        peak_freq = peak_indices * Fs / window_size;
% Multiply peak indices by 25
convertindex_peak_freq = peak_indices * 25;
sorted_convertindex_peak_freq = sort(convertindex_peak_freq);
% Find the characteracter with matching frequencies
        same_character = find_same_character(sorted_convertindex_peak_freq, frq_characteracter);
    Text_T.String = ['Matched characteracter: ', same_character];
        % Store the decoded characteracter
        decoded_characteracters{frame_index} = same_character;
    end
    % Concatenate the decoded characteracters
    decoded_string = strjoin(decoded_characteracters, '');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function same_character = find_same_character(frequencies, frq_characteracter)
    same_character = '';  
    for character = fieldnames(frq_characteracter)' % Check each characteracter
        % Compare frequencies
        if (frq_characteracter.(character{1})(1) == frequencies(1)) && ...
           (frq_characteracter.(character{1})(2) == frequencies(2)) && ...
           (frq_characteracter.(character{1})(3) == frequencies(3))
            if character{1} == 'space' 
                           same_character = ' '; 
            else
            same_character = character{1};
            break; % No need to check further once a match is found
            end
        end
    end
    if isempty(same_character)
        same_character = 'y'; % Assign a placeholder if no match found
    end
        Text_T.String = ['same character: ', same_character];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to get default characteracter frequencies
function frq_characteracter = getDefaultfrq_characteracter()
   frq_characteracter.a = [100, 1100, 2500]; 
frq_characteracter.b = [100, 1100, 3000];
frq_characteracter.c = [100, 1100, 3500];
frq_characteracter.d = [100, 1300, 2500];
frq_characteracter.e = [100, 1300, 3000];
frq_characteracter.f = [100, 1300, 3500];
frq_characteracter.g = [100, 1500, 2500];
frq_characteracter.h = [100, 1500, 3000];
frq_characteracter.i = [100, 1500, 3500];
frq_characteracter.j = [300, 1100, 2500];
frq_characteracter.k = [300, 1100, 3000];
frq_characteracter.l = [300, 1100, 3500];
frq_characteracter.m = [300, 1300, 2500];
frq_characteracter.n = [300, 1300, 3000];
frq_characteracter.o = [300, 1300, 3500];
frq_characteracter.p = [300, 1500, 2500];
frq_characteracter.q = [300, 1500, 3000];
frq_characteracter.r = [300, 1500, 3500];
frq_characteracter.s = [500, 1100, 2500];
frq_characteracter.t = [500, 1100, 3000];
frq_characteracter.u = [500, 1100, 3500];
frq_characteracter.v = [500, 1300, 2500];
frq_characteracter.w = [500, 1300, 3000];
frq_characteracter.x = [500, 1300, 3500];
frq_characteracter.y = [500, 1500, 2500];
frq_characteracter.z = [500, 1500, 3000];
frq_characteracter.space = [500, 1500, 3500];
end

