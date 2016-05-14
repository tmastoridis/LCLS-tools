function sys_prm = model_init()

sys_prm = struct();

sys_prm.frf = 1.3e9;
sys_prm.wrf = 2*pi*sys_prm.frf;

sys_prm.cavity = struct();
sys_prm.cavity.r_Q = 1036;
sys_prm.cavity.dw_true = -1e3 * 2 * pi;
sys_prm.cavity.QL_true = 4.12e7;
sys_prm.cavity.Gn_true = 1.0;