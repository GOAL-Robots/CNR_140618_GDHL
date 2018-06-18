function y = dwsnl(time_d,decay_rate,coeff_K)

y = -exp(time_d/decay_rate(2)).*coeff_K^2.*decay_rate(1).*(decay_rate(2)*(-decay_rate(1)+decay_rate(2))+time_d*(decay_rate(1)+decay_rate(2)))/(decay_rate(1)+decay_rate(2))^3;
% dwsnh exp(-((time_d+decay_rate(1)+decay_rate(2))/decay_rate(1)))*coeff_K^2*decay_rate(1)*(time_d*(decay_rate(1)+decay_rate(2))+decay_rate(2)*(3*decay_rate(1)+decay_rate(2)))/(decay_rate(1)+decay_rate(2))^3;
% dwpsl -exp(-((-time_d+decay_rate(1)+decay_rate(2))/decay_rate(2)))*coeff_K^2*(exp(1+decay_rate(1)/decay_rate(2))*decay_rate(1)*(decay_rate(2)*(-decay_rate(1)+decay_rate(2))+time_d*(decay_rate(1)+decay_rate(2)))+decay_rate(2)*(time_d*(decay_rate(1)+decay_rate(2))-decay_rate(1)*(decay_rate(1)+3*decay_rate(2))))/(decay_rate(1)+decay_rate(2))^3;
% dwpsh exp(-1-time_d/decay_rate(1)-decay_rate(1)/decay_rate(2))*coeff_K^2*decay_rate(2)*(exp(1+decay_rate(1)/decay_rate(2))*(decay_rate(1)*(decay_rate(1)-decay_rate(2))-time_d*(decay_rate(1)+decay_rate(2)))+exp(time_d*(1/decay_rate(1)+1/decay_rate(2)))*(-time_d*(decay_rate(1)+decay_rate(2))+decay_rate(1)*(decay_rate(1)+3*decay_rate(2))))/(decay_rate(1)+decay_rate(2))^3;
% dwnsl exp(-(-time_d+decay_rate(1)+decay_rate(2))/decay_rate(2))*coeff_K^2*decay_rate(2)*(-time_d*(decay_rate(1)+decay_rate(2))+decay_rate(1)*(decay_rate(1)+3*decay_rate(2)))/(decay_rate(1)+decay_rate(2))^3;
% dwnsh exp(-time_d/decay_rate(1))*coeff_K^2*decay_rate(2)*(decay_rate(1)*(-decay_rate(1)+decay_rate(2))+time_d*(decay_rate(1)+decay_rate(2)))/(decay_rate(1)+decay_rate(2))^3;
