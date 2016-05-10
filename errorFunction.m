function er = errorFunction(x, sys_prm, data, w)

Hc = cavity_model(sys_prm, x(1), x(2), w);

er = sum(abs(Hc-data).^2);

