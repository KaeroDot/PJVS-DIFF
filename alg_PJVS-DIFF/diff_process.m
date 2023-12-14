% Inputs:
% sigconfig - structure with configuration data:
%   .fs - sampling frequency (Hz)
%   .fseg - frequency of changing PJVS steps (Hz)
%   .PRs - masked after PJVS step change (samples)
%   .PRe - masked before PJVS step change (samples)
%   .MRs - masked after record start (samples)
%   .MRe - masked before record start (samples)
% y - sampled data
% Uref1period - reference values of PJVS voltages for one PJVS period
% dbg - debug structure or debug level
%
% Outputs:
% A_rms_total - RMS amplitude of the DUT signal calculated from whole data
% A_rms_t - RMS amplitude calculated from every signal period
% A_fft - amplitude of the DUT signal calculated using FFT (no special window)
% A_fft_l - noise amplitude to the left of the DUT signal calculated using FFT (no special window)
% A_fft_r - noise amplitude to the right amplitude of the DUT signal calculated using FFT (no special window)

function [A_rms_total, A_rms_t, A_fft, A_fft_l, A_fft_r] = diff_process(sigconfig, y, Uref1period, dbg) %<<<1

    % Check inputs %<<<2
    if nargin ~= 4
        error('diff_process: bad number of input arguments!')
    end
    dbg = check_gen_dbg(dbg);
    sigconfig = check_sigconfig(sigconfig);
    % ensure the directory for plots exists
    if dbg.v && (dbg.saveplotsfig || dbg.saveplotspng)
        if ~exist(dbg.plotpath, 'dir')
            mkdir(dbg.plotpath);
        end
    end

    % Calc and check properties %<<<2
    % Length of PJVS segments in samples (samples between two different PJVS steps):
    segmentlen = sigconfig.fs./sigconfig.fseg;
    if rem(segmentlen, 1)
        warning('diff_process: sampling frequency and PJVS segments frequency is not divisible, rounding errors can occur!')
        segmentlen = round(segmentlen);
    end

    % Number of segments in a period (PJVS steps in a signal periods):
    segments = sigconfig.fseg./sigconfig.f;
    if rem(segments, 1)
        warning('diff_process: number of PJVS segments in a signal period is not divisible, rounding errors can occur!')
        segments = round(segments);
    end

    % Check that at least one segment will be left after MRs MRe:
    if (sigconfig.MRs + sigconfig.MRe >= numel(y) + segmentlen)
        error(sprintf('diff_process: no segments left in data after masking MRs, MRe: MRs(%g) + MRe(%g) > segmentlen(%g) + size of y(%g).',...
            sigconfig.MRs,...
            sigconfig.MRe,...
            segmentlen,...
            numel(y)))
    end % if

    % Check that something will be left after PRs PRe:
    if (sigconfig.PRs + sigconfig.PRe >= segmentlen)
        error(sprintf('diff_process: no samples left in segments after masking PRs, PRe: PRs(%g) + PRe(%g) >= segmentlen(%g).',...
            sigconfig.PRs,...
            sigconfig.PRe,...
            segmentlen))
    end % if

    % Identify PJVS segments %<<<2
    % Find out indexes with start of PJVS segments in the sampled signal:
    % set section values because the function is copied from QPSW and needs it:
    dbg.section = [0 0];
    Spjvs = pjvs_ident_segments(y, sigconfig.MRs, sigconfig.MRe, segmentlen, dbg);

    % If first segment is not whole (i.e. samples in the first segment is not equal
    % to segmentlen), set sigconfig.MRs so this part will be removed by function
    % pjvs_split_segments. 
    if Spjvs(2) - Spjvs(1) ~= Spjvs(3) - Spjvs(2)
        sigconfig.MRs = Spjvs(2);
    end
    % If last segment is not whole (i.e. samples in the last segment is not equal
    % to segmentlen), set sigconfig.MRe.
    if Spjvs(end) - Spjvs(end-1) ~= Spjvs(end-1) - Spjvs(end-2)
        sigconfig.MRe = Spjvs(end) - Spjvs(end-1) - 1;
    end

    % Split signal into PJVS segments %<<<2
    [s_y, s_mean, s_std, s_uA, s_slope, neglected_segments] = pjvs_split_segments(y, Spjvs, sigconfig.MRs, sigconfig.MRe, sigconfig.PRs, sigconfig.PRe, dbg, 1);
    % If number of segments is not integer multiple of signal periods, remove
    % last extra segments:
    if rem(size(s_y, 2), segments) % XXX divide by segment len? no, divide by number of segments in period!
        disp('Removing extra segments to keep integer multiple of signal periods.')
        cut = floor(size(s_y, 2)/segments).*segments;
        s_y = s_y(:, 1:cut);
        s_mean = s_mean(:, 1:cut);
        s_std = s_std(:, 1:cut);
        s_uA = s_uA(:, 1:cut);
        s_slope = s_slope(:, 1:cut);
    end % if ~isempty(neglected_segments)

    % Identify phase of the signal and match PJVS reference values to all segments in whole signal %<<<2
    % Because in differential measurement every segment got similar s_mean, it
    % is not used to identify phase of PJVS vs DUT signal. Instead slope of
    % segment is used.
    Uref = pjvs_ident_Uref(s_slope, Uref1period, dbg);
    % add 90 deg phase because smallest slope is top of waveform, not minimum, and pjvs_ident_uref calculates to minimum
    % (XXX this is very bad method)
    Uref = circshift(Uref, round(segments./4));

    plot_sine_fit(s_y, Uref, segments, sigconfig.fs, sigconfig.f, dbg);

    % Calculate RMS %<<<2
    % Calculate total RMS amplitude for the whole signal
    % Because averages per PJVS steps are not calculated, sinc correction is not
    % needed, and RMS can be calculated directly. The s_y should already contain
    % only integer multiple of signal periods.
    % (rms function exist, but in matlab it is part of some costly package)
    allsegments = (s_y + Uref)(:);
    A_rms_total = mean(allsegments.^2).^0.5;

    % Calculate RMS values for every single signal period:
    N = 1; % RMS every 1 period
    % number of samples in a single period
    samplesinNperiods = N.*segments.*segmentlen;
    % just to be sure that integer number of periods is used:
    integercut = N.*segments.*floor(size(s_y, 2)./(N.*segments));
    y_real = [s_y(:, 1:integercut) + Uref(1:integercut)];
    % reshape samples by one period:
    tmp = reshape(y_real(:), samplesinNperiods, []);
    A_rms_t = mean(tmp.^2, 1).^0.5;
    % plot RMS values
    plot_rms(A_rms_t, A_rms_total, dbg);

    % Calculate FFT %<<<2
    % Calculate FFT of whole signal.
    % frequency axis is not correct. Due to PRs, PRe, the output of SP-WFFT
    % algorithm has to be corrected by some factor.
    DI.fs.v = 1;
    DI.y.v = y_real(:)';
    DI.y.v(isnan(DI.y.v)) = [];
    DOfft = qwtb('SP-WFFT', DI);
    % Find signal amplitude as highest peak:
    [A_fft, ind] = max(DOfft.A.v);
    A_fft = A_fft(1);
    f_fft = DOfft.f.v(ind(1));
    A_fft_l = DOfft.A.v(ind(1)-1);
    A_fft_r = DOfft.A.v(ind(1)+1);
    % plot spectrum
    plot_fft(DOfft, f_fft, A_fft, A_fft_l, A_fft_r, dbg);

    % Calculate amplitudes using fitting  %<<<2
    % % This gives quite bad results compared to simple RMS.
    % % Even for ideal signal and 100 periods the error is in few ppm.
    % % Calculate amplitude of signal using FPNLSF
    % % Amlitude will be calculated for every N signal periods:
    % N = 100;
    % N_vec = [N.*segments : N.*segments : size(s_y, 2)];
    % N_vec = [1 N_vec];
    % if N_vec(end) - N_vec(end-1) ~= N.*segments;
    %     N_vec(end) = [];
    % end % if
    %
    % t = [0 : prod(size(s_y))-1]./sigconfig.fs;
    % t = reshape(t, size(s_y));
    % DI.fs.v = sigconfig.fs;
    % DI.fest.v = sigconfig.f;
    % for j = 1 : numel(N_vec)-1
    %     DI.y.v = s_y(:, N_vec(j) : N_vec(j+1) - 1) + Uref(N_vec(j) : N_vec(j+1) - 1);
    %     DI.y.v = DI.y.v(:);
    %     DI.t.v = t(:, N_vec(j) : N_vec(j+1) - 1);
    %     DI.t.v = DI.t.v(:);
    %     FPNLSF(j) = qwtb('FPNLSF', DI);
    %     RMS(j) = mean(DI.y.v.^2).^0.5;
    % end % for j = 1 : numel(N_vec)
    % amps = [[FPNLSF.A].v];
    % figure
    % hold on
    % plot([amps./sqrt(2) - A_rms_total].*1e6, '-k')
    % hold off

end % function diff_process

function plot_fft(DOfft, f_fft, A_fft, A_fft_l, A_fft_r, dbg) %<<<1
    % Plot spectrum
    if dbg.v
        if dbg.rms_plot
            figure('visible',dbg.showplots)
            hold on
            loglog(DOfft.f.v, DOfft.A.v)
            indA = find(DOfft.A.v == A_fft);
            plot(DOfft.f.v(indA(1)), DOfft.A.v(indA(1)), 'or')
            indl = find(DOfft.A.v == A_fft_l);
            plot(DOfft.f.v(indl(1)), DOfft.A.v(indl(1)), 'ok')
            indr = find(DOfft.A.v == A_fft_r);
            plot(DOfft.f.v(indr(1)), DOfft.A.v(indr(1)), 'ok')
            xlabel('f (a.u.)')
            ylabel('U (V)')
            title(sprintf(...
                'FFT spectrum, unscaled x-axis\nnoise: left %.3g dB, right %.3g dB', ...
                20*log10(DOfft.A.v(indl(1))/DOfft.A.v(indA(1))), ...
                20*log10(DOfft.A.v(indr(1))/DOfft.A.v(indA(1)))));
            hold off
            fn = fullfile(dbg.plotpath, 'fft_spectrum');
            if dbg.saveplotsfig saveas(gcf(), [fn '.fig'], 'fig') end
            if dbg.saveplotspng saveas(gcf(), [fn '.png'], 'png') end
            close
        end % if dbg.rms_plot
    end % if dbg.v
end % function plot_rms

function plot_sine_fit(s_y, Uref, segments, fs, fnom, dbg) %<<<1
    % Plot segments for first and last period. First plot shows segments with
    % sine fit, second plot shows segments with sine subtracted.
    if dbg.sine_fit
        % s_y indexes of first period and last period:
        ind = {[1, segments],...
               [size(s_y, 2) - segments + 1, size(s_y, 2)]};
        fn_append = {'first', 'last'};
        for j = 1:numel(ind)
            % s_y indexes for actual period
            id = ind{j};
            % cut one period of signal
            % selected segments:
            selsegments = (s_y(:, id(1):id(2)) + Uref(id(1):id(2)));
            % segments contains all samples, while bad samples (sigconfig.PRs/e and sigconfig.MRs/e) are
            % represented as NaN values.
            % absolute time vector:
                                            % t = [0:numel(selsegments(:))-1]./fs;
            t_offset = (id(1) - 1).*numel(s_y(:,id(1)))./fs;
            t = [0:numel(selsegments(:))-1]./fs;
            % remove nan values for fitting and prepare for FPNLSF:
            % FPNLSF have to be used because this algorithm can handle many
            % discontinuities in the signal (but cannot handle NaN values).
            DI.y.v = selsegments(:);
            % Fitting does not like starting from nonzero values.
            % When you fit the same data with t(1) = 0, it gets much better results than if t(1) > 0!
            % Therefore for fit, t_offset is removed, and later the phase will
            % be recalculated to make proper plot data
            DI.t.v = t;
            idx = isnan(DI.y.v);
            DI.t.v(idx) = [];
            DI.y.v(idx) = [];
            DI.fs.v = fs;
            DI.fest.v = fnom;
            Fit = qwtb('FPNLSF', DI);

            % reference plot of the fit:
            % reconstruct full sine waveform
            ysin = Fit.A.v.*sin(2*pi*Fit.f.v*t + Fit.ph.v);
            figure('visible',dbg.showplots)
            hold on
            plot(DI.t.v + t_offset, DI.y.v, 'x', t + t_offset, ysin, '-');
            xlabel('t (s)')
            ylabel('U (V)')
            title(sprintf('Sampled data with PJVS reference value fitted by sine wave\n%s period', fn_append{j}))
            legend('samples', 'fit by sine', 'location','southoutside','orientation', 'horizontal')
            hold off
            fn = fullfile(dbg.plotpath, ['signal_sine_fit_period_' fn_append{j}]);
            if dbg.saveplotsfig saveas(gcf(), [fn '.fig'], 'fig') end
            if dbg.saveplotspng saveas(gcf(), [fn '.png'], 'png') end
            close

            % reference plot of the segments:
            % reshape sine waveform for subtracting from segments:
            ysinseg = reshape(ysin, size(selsegments));
            figure('visible',dbg.showplots)
            % voltage offset for plotting many line plots:
            offset = 3.*mean(std(selsegments - ysinseg, 1));
            % round offset to tenths of uV:
            offset = 1e-5.*ceil(1e5.*offset);
            leg = {};
            hold on
            for k = 1:size(selsegments, 2)
                % plot with increasing voltage offset
                plot(selsegments(:,k) - ysinseg(:,k) + (k-1).*offset);
                leg{end+1} = num2str(k);
            end % for
            xlabel('samples')
            ylabel('U (V)')
            title(sprintf('PJVS selected segments minus fit sine wave,\n offseted by %g uV', 1e6.*offset))
            legend(leg, 'location', 'eastoutside');
            hold off
            fn = fullfile(dbg.plotpath, ['signal_minus_fit_period_' fn_append{j}]);
            if dbg.saveplotsfig saveas(gcf(), [fn '.fig'], 'fig') end
            if dbg.saveplotspng saveas(gcf(), [fn '.png'], 'png') end
            close
        end % for
    end % if
end % function plot_sine_fit(s_y, Uref, segments, fs, dbg)

function plot_rms(A_rms_t, A_rms_total, dbg) %<<<1
% plot RMS values
    if dbg.rms_values
        figure('visible', dbg.showplots)
        hold on
        plot([A_rms_t - A_rms_total].*1e6, '-');
        xlabel('signal period')
        ylabel('U (uV), RMS')
        title(sprintf('RMS values minus total RMS of the whole signal\nRMS_{total} = %.9f V\n std(RSM) = %.3f uV', A_rms_total, 1e6.*std(A_rms_t)))
        hold off
        fn = fullfile(dbg.plotpath, ['rms_values_per_period']);
        if dbg.saveplotsfig saveas(gcf(), [fn '.fig'], 'fig') end
        if dbg.saveplotspng saveas(gcf(), [fn '.png'], 'png') end
        close
    end % if dbg.v
end % function plot_rms

% vim settings modeline: vim: foldmarker=%<<<,%>>> fdm=marker fen ft=octave textwidth=80 tabstop=4 shiftwidth=4
