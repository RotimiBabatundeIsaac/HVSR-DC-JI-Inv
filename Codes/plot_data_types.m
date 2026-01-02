% Load and verify
load('DATA/sta_datapack.mat');

% Check a few stations
for i = [1, 45, 91]
    fprintf('\n=== %s (%.3f, %.3f) ===\n', sta(i).name, sta(i).lat, sta(i).lon);
    fprintf('  GRV: %d periods, U = [%.2f - %.2f] km/s\n', ...
        numel(sta(i).ugr.T), min(sta(i).ugr.obs), max(sta(i).ugr.obs));
    fprintf('  PHV: %d periods, c = [%.2f - %.2f] km/s\n', ...
        numel(sta(i).cph.T), min(sta(i).cph.obs), max(sta(i).cph.obs));
    fprintf('  Ellip: %d periods\n', numel(sta(i).ellip.T));
    fprintf('  HVSR: %d freqs\n', numel(sta(i).hvsr.f));
end

% Quick plot for one station
figure;
idx = 45;  % Pick a station with all data types
subplot(2,2,1); plot(sta(idx).ugr.T, sta(idx).ugr.obs, 'o-'); 
xlabel('Period (s)'); ylabel('U (km/s)'); title('Group Velocity');

subplot(2,2,2); plot(sta(idx).cph.T, sta(idx).cph.obs, 's-');
xlabel('Period (s)'); ylabel('c (km/s)'); title('Phase Velocity');

subplot(2,2,3); plot(sta(idx).ellip.T, sta(idx).ellip.obs, '^-');
xlabel('Period (s)'); ylabel('H/V'); title('Ellipticity');

subplot(2,2,4); semilogx(sta(idx).hvsr.f, sta(idx).hvsr.obs, '-');
xlabel('Frequency (Hz)'); ylabel('H/V'); title('HVSR');

sgtitle(sta(idx).name);