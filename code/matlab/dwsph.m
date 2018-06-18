function y = dwsph(time_d,decay_rate,coeff_K)

y = exp(-(time_d+decay_rate(1)+decay_rate(2))/decay_rate(1)).*coeff_K^2.*(exp(1+decay_rate(2)/decay_rate(1))*decay_rate(2)*(decay_rate(1)*(-decay_rate(1)+decay_rate(2))+time_d*(decay_rate(1)+decay_rate(2)))+decay_rate(1)*(time_d*(decay_rate(1)+decay_rate(2))+decay_rate(2)*(3*decay_rate(1)+decay_rate(2))))/(decay_rate(1)+decay_rate(2))^3;