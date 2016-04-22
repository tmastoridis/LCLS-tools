function sys_prm = model_init()

sys_prm = struct();

sys_prm.frf = 1.3e9;
sys_prm.wrf = 2*pi*sys_prm.frf;
