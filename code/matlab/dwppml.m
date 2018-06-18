function y = dwppml(time_d,decay_rate,coeff_K)

y = exp(-(-time_d+decay_rate(1)+decay_rate(2))/decay_rate(2))*coeff_K^2*(decay_rate(2)*(-decay_rate(1)^2-2*decay_rate(1)*decay_rate(2)+decay_rate(2)^2 +time_d*(decay_rate(1)+decay_rate(2)))+exp(1+decay_rate(1)/decay_rate(2))*decay_rate(1)*(2*decay_rate(2)^2+time_d*(decay_rate(1)+decay_rate(2))))/(decay_rate(2)*(decay_rate(1)+decay_rate(2))^3);
