function y = dwppmr(time_d,decay_rate,coeff_K)

y = exp(-(time_d+decay_rate(1)+decay_rate(2))/decay_rate(1))*coeff_K^2*(-time_d*(decay_rate(1)+decay_rate(2))*(decay_rate(1)+exp(1+decay_rate(2)/decay_rate(1))*decay_rate(2))+decay_rate(1)*(decay_rate(1)^2+2*(-1+exp(1+decay_rate(2)/decay_rate(1)))*decay_rate(1)*decay_rate(2)-decay_rate(2)^2))/(decay_rate(1)*(decay_rate(1)+decay_rate(2))^3);