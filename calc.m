addpath('~/metrologie/Q-Wave/qwtb/qwtb')
addpath('alg_PJVS-DIFF/')

% % -------------- simulation
% segments = 20; % XXX fity nefunguji pro jine nez 20!
% fs_error = 0;
%
% % test generatoru:
% sigconfig.f = 60;
% sigconfig.A = 7;
% sigconfig.ph = 0;
% sigconfig.fs = 6*1e6;
% sigconfig.Lm = 10;
% sigconfig.noise = 2e-6;
% sigconfig.fseg = segments .* sigconfig.f;
% sigconfig.fm = 75e9;
% sigconfig.apply_filter = 1;
%
% dbg = 0;
% dbg = check_gen_dbg(dbg);
% dbg.section = [0 0];
% sigconfig.fs = sigconfig.fs + fs_error;
% sigconfig.MRs = 0;
% sigconfig.MRe = 0;
% sigconfig.PRs = 200;
% sigconfig.PRe = 200;
%
%         % [D, Uref, Uref1period, Sid] = diff_simulator(sigconfig, dbg);
%         % [RMStotal, RMSvec, A_fft, Al_fft, Ar_fft] = diff_process(sigconfig, D, Uref1period, dbg);
%         %
%         % (RMStotal - sigconfig.f.*2^-0.5)./(sigconfig.f.*2^-0.5)
%         % (A_fft - sigconfig.f.*2^-0.5)./(sigconfig.f.*2^-0.5)
%
% % qwtb alg method:
% DI.fs.v = sigconfig.fs;
% DI.fseg.v = sigconfig.fseg;
% DI.f.v = sigconfig.f;
% DI.y.v = D;
% DI.Uref.v = Uref1period;
% DI.Rs.v = sigconfig.PRs;
% DI.Re.v = sigconfig.PRe;
% DI.Ms.v = sigconfig.MRs;
% DI.Me.v = sigconfig.MRe;
% DI.plots.v = 5;
% DI.data_folder.v = 'QPSW_plots';
%
% DO = qwtb('PJVS-DIFF', DI);
%
% (DO.U.v - sigconfig.f.*2^-0.5)./(sigconfig.f.*2^-0.5)

% -------------- real data calculation 2

DI.fs.v = 3200*1e3; % urceno z poctu vzorku v 1 periode v datech
DI.fseg.v = 20.*DI.f.v; % urceno z poctu segmentu v 1 periode v datech
% DI.y.v = load('mereni/data copy z ACVScontrol jedna perioda/f5720A_1V_rms_1kHz_1diffperiod_from_ACSVcontrol_001.csv');
% DI.y.v = load('mereni/data copy z ACVScontrol jedna perioda/f5720A_1V_rms_1kHz_1diffperiod_from_ACSVcontrol_002.csv');
DI.y.v = load('mereni/data copy z ACVScontrol jedna perioda/f5720A_1V_rms_1kHz_1diffperiod_from_ACSVcontrol_003.csv');
% doufam ze to je spravne:
DI.Uref.v = [0 0.437021175 0.831205621 1.14404355 1.34495791 1.41418898 1.34495791 1.14404355 0.831205621 0.437021175 0 -0.437021175 -0.831205621 -1.14404355 -1.34495791 -1.41418898 -1.34495791 -1.14404355 -0.831205621 -0.437021175];
DI.Rs.v = 0;
DI.Re.v = 0;
DI.Ms.v = 0;
DI.Me.v = 0;
DI.plots.v = 5;
DI.data_folder.v = 'QPSW_plots';

DO = qwtb('PJVS-DIFF', DI);

(DO.U.v - 2^-0.5)./(2^-0.5)

% % % -------------- real data calculation 1
% % "diferencialni mereni AC s PJVS exetrnal control
% % vlastni NI5622
% % 10 MHz z PJVS
% % zem racku na PJVS
% % trigger level z PJVS, ale to nefunguje nevim proc
% % kalibrator F5720A 9385202
% % open guard, open lo,
% % phase lock z PJVS
% % 1 V 1 kHz
% % PJVS hodnoty podle external control:
% % 0
% % 0.437021175
% % 0.831205621
% % 1.14404355
% % 1.34495791
% % 1.41418898
% % 1.34495791
% % 1.14404355
% % 0.831205621
% % 0.437021175
% % 0
% % -0.437021175
% % -0.831205621
% % -1.14404355
% % -1.34495791
% % -1.41418898
% % -1.34495791
% % -1.14404355
% % -0.831205621
% % -0.437021175
% % "
%
% % samples count:: 6000000
% % sampling rate [Sa/s]:: 6000000.00000000
% %startmatrix:: record sample data gains [V]
% data_gain = 5.5250876E-9;
%
% Uref1period = [0, ...
%     0.437021175, ...
%     0.831205621, ...
%     1.14404355, ...
%     1.34495791, ...
%     1.41418898, ...
%     1.34495791, ...
%     1.14404355, ...
%     0.831205621, ...
%     0.437021175, ...
%     0, ...
%     -0.437021175, ...
%     -0.831205621, ...
%     -1.14404355, ...
%     -1.34495791, ...
%     -1.41418898, ...
%     -1.34495791, ...
%     -1.14404355, ...
%     -0.831205621, ...
%     -0.437021175];
%
% fnom = 1e3; % 1 kHz signal, nominal value
% sigconfig.fs = 6e6; % sampling frequency
% segments = 20;
% sigconfig.fseg = fnom.*segments; % sampling frequency
% % segmentlen = fs/(fnom.*20);
%
% % load('../manual/230728_AC_differential_10s_0002/RAW/G0001-A0001.mat')
% load('../manual/230728_AC_differential_1s_0001/RAW/G0001-A0001.mat')
% % load('../manual/230728_AC_differential_1s_0004/RAW/G0001-A0001.mat')
% % load('../manual/230728_AC_differential_333.333s_0003/RAW/G0001-A0001.mat')
% % load('G0001-A0001.mat')
% y = y.*data_gain;
% sigconfig.MRs = 0; % this spoils diff. calculations if together not integer number of samples in period!
% sigconfig.MRe = 0; % this spoils diff. calculations if together not integer number of samples in period!
% sigconfig.PRs = 80;
% sigconfig.PRe = 40;
%
% dbg = check_gen_dbg(5);
% dbg.section = [0 0];
%
% [RMStotal, RMSvec] = diff_process(sigconfig, y, [], [], Uref1period, [], [], dbg);
%
%
