# ======================================================================
# PURPOSE:
# Read and process THEMIS FGM data (CDF) to identify particle injections
# using magnetic field signatures (dipolarization).
# ======================================================================

# %% ---------------------------- IMPORTS ----------------------------
import datetime
import os
import cdflib
from cdflib import cdfepoch
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


# %% ---------------------------- 1. TIME INTERVAL -------------------
start_time = datetime.datetime(2007, 3, 23, 0, 0)
end_time   = datetime.datetime(2007, 3, 24, 0, 0)

start_plot_dt = datetime.datetime(2007, 3, 23, 11, 0)
end_plot_dt   = datetime.datetime(2007, 3, 23, 12, 0)


# %% ---------------------------- 2. DIRECTORY -----------------------
#FGM instrument
directory_FGM = r'C:\Users\karen\Programas\particle_injection_projects\Data\FGM'
thed_file_FGM = 'thd_l2_fgm_20070323_v01.cdf'
path_thed_FGM = os.path.join(directory_FGM, thed_file_FGM)

#SST instrument
directory_SST = r'C:\Users\karen\Programas\particle_injection_projects\Data\SST'
thed_file_SST = 'thd_l2_sst_20070323_v01.cdf'
path_thed_SST = os.path.join(directory_SST, thed_file_SST)


# %% ---------------------------- 3. OPEN CDF ------------------------
#FGM instrument
cdf_thed_FGM = cdflib.CDF(path_thed_FGM)
print("Variables in CDF:", cdf_thed_FGM.cdf_info().zVariables)

#SST instrument
cdf_thed_SST = cdflib.CDF(path_thed_SST)
print("Variables in CDF:", cdf_thed_SST.cdf_info().zVariables)


# %% ---------------------------- 4. READ VARIABLES FGM ------------------
#FGM instrument
epoch_raw_FGM = cdf_thed_FGM.varget('thd_fgs_time')  # em segundos desde 1970
Btotal    = cdf_thed_FGM.varget('thd_fgs_btotal')
Bgse      = cdf_thed_FGM.varget('thd_fgs_gse')

# sanity check
print("Epoch sample:", epoch_raw_FGM[:5])
print("Btotal sample:", Btotal[:5])
print("Bgse sample:\n", Bgse[:5, :])

# Covert epoch to datetime
# Unix epoch (seconds since 1970-01-01)
epoch_dt_FGM = pd.to_datetime(epoch_raw_FGM, unit='s', origin='unix')

# create DataFrame
df_thed_FGM = pd.DataFrame(
    {
        'Btotal': Btotal,
        'Bx': Bgse[:, 0],
        'By': Bgse[:, 1],
        'Bz': Bgse[:, 2],
    },
    index=pd.DatetimeIndex(epoch_dt_FGM)
)
df_thed_FGM.index.name = 'time'

# check
print(df_thed_FGM.index.min())
print(df_thed_FGM.index.max())


# %% ---------------------------- 5. READ VARIABLES SST ------------------
# SST instrument (PSIF)
epoch0 = cdf_thed_SST.varget('thd_psif_epoch0')           # base CDF epoch
delta_time = cdf_thed_SST.varget('thd_psif_delta_time')   # segundos desde epoch0
spec_ion_flux = cdf_thed_SST.varget('thd_psif_en_eflux')  # fluxos espectrais series_ion_flux = cdf_thed_SST.varget('thd_psif_flux')    # fluxos séries
series_ion_flux = cdf_thed_SST.varget('thd_psif_flux')    # fluxos séries


# -------------------- INSPECTION --------------------

variables = {
    'epoch0': epoch0,
    'delta_time': delta_time,
    'spec_ion_flux': spec_ion_flux,
    'series_ion_flux': series_ion_flux
}

for name, var in variables.items():
    print(f"\n{name.upper()}:")
    print(f"  Tipo: {type(var)}")
    try:
        print(f"  Shape: {np.shape(var)}")
    except:
        print("  Shape: não aplicável")
    
    # imprimir os primeiros elementos
    if isinstance(var, (np.ndarray, list)):
        print("  Primeiros valores:", np.array(var).flatten()[:5])
    else:
        print("  Valor:", var)

#%% -------------------- CONVERT EPOCH TO DATETIME --------------------
# eixo temporal completo do SST
epoch0_dt = pd.Timestamp(cdfepoch.to_datetime(epoch0))
delta_td = pd.to_timedelta(delta_time, unit='s')

time_sst = epoch0_dt + delta_td

print(time_sst[0], time_sst[-1])

# ESTÁ DANDO ERRO AQUI, PRECISO VERIFICAR DEPOIS

# %% ---------------------------- 6. FILTER TIME INTERVAL -----------
df_plot = df_thed_FGM[
    (df_thed_FGM.index >= start_plot_dt) &
    (df_thed_FGM.index <= end_plot_dt)
]
epoch = df_plot.index

# %% ---------------------------- 7. PLOTS ---------------------------
fig, axs = plt.subplots(2, 1, figsize=(12, 6), sharex=True)
fig.subplots_adjust(hspace=0.05)

# (a) |B|
axs[0].plot(epoch, df_plot['Btotal'], color='black')
axs[0].set_ylabel('|B| (nT)', fontsize=13)
axs[0].grid(True)
axs[0].tick_params(axis='both', labelsize=11)
axs[0].text(0.02, 0.85, '(a)', transform=axs[0].transAxes,
            fontsize=13, fontweight='bold')

# (b) B components (GSE)
axs[1].plot(epoch, df_plot['Bx'], label='Bx (GSE)')
axs[1].plot(epoch, df_plot['By'], label='By (GSE)')
axs[1].plot(epoch, df_plot['Bz'], label='Bz (GSE)')
axs[1].set_ylabel('B (nT)', fontsize=13)
axs[1].set_xlabel('Time (UT)', fontsize=13)
axs[1].grid(True)
axs[1].legend(fontsize=11)
axs[1].tick_params(axis='both', labelsize=11)
axs[1].text(0.02, 0.85, '(b)', transform=axs[1].transAxes,
            fontsize=13, fontweight='bold')

plt.show()

# %%
