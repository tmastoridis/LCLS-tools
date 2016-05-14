function er = cavity_error(sys_prm, x, w, data)

Hc = cavity_model(sys_prm, x(1), x(2), x(3), w);

er = sum(abs(Hc-data).^2);

