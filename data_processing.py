import tqdm
import numpy as np
import pandas as pd
import warnings
warnings.filterwarnings("ignore")

##################################################
# Prepare lab values and self-defined endpoints
##################################################

# link to data dict
# https://docs.google.com/spreadsheets/d/1qpa9KFp36x1qQff14OEhtviQqrJQO-qQIKuqt-QX_qA/edit#gid=0
# ep_path = '/data/processed_data/detailed_longitudinal/R10/backup/splits_2024_04_01/'
ep_path = '/data/processed_data/detailed_longitudinal/R10/detailed_longitudinal_2024-01-24.csv'
# Andrius file
# '/data/processed_data/detailed_longitudinal/R10/detailed_longitudinal_DF10_2022-11-11.csv'
lab_path = '/data/processed_data/kela_lab/kanta_lab_v2_2023-12-12.csv.gz'

ep_cols = ['FINREGISTRYID','EVENT_AGE','EVENT_DAY','ICDVER', 'E102','N083','2503B']
lab_cols = ['FINREGISTRYID', 'LAB_DATE_TIME', 'LAB_VALUE', 'LAB_UNIT', 'OMOP_ID', 'LAB_ABNORMALITY',
            'REFERENCE_VALUE_TEXT']

# ICD10: E102&N083|E102
# ICD9: 2503B

reader = pd.read_csv(ep_path, iterator=True)#, chunksize=10**6)
ep_df = pd.DataFrame(columns=ep_cols)

N = 2487878501
n = 0

while n < N:
    tmp = reader.get_chunk(500000)
    n += 500000
    tmp['E102'] = (tmp == 'E102').any(axis=1)
    tmp['N083'] = (tmp == 'N083').any(axis=1)
    tmp['2503B'] = (tmp == '2503B').any(axis=1)
    tmp = tmp[(tmp[['E102','N083','2503B']] == True).any(axis=1) == True]
    tmp = tmp[ep_cols]
    ep_df = pd.concat([ep_df,tmp])
    if n in [1000000, 10000000, 50000000, 100000000, 500000000,
             1000000000, 1500000000, 2000000000, 2487000000]:
        ep_df.to_csv('df_nephropathy.csv', index=None)
        print(n)
ep_df.to_csv('df_nephropathy.csv', index=None)


# 3557 - 3050449: Albumin [Mass/time] in Urine collected for unspecified duration
    # unit: ¬µg/min
# 4836 - 3049506: Microalbumin [Mass/time] in Urine collected for unspecified duration
    # unit: ¬µg/min
# 2513 - 3020876: Protein [Mass/time] in 24 hour Urine
    # unit: mg
# 3020564: Creatinine [Moles/volume] in Serum or Plasma
    # unit: ¬µmol or ¬µmol/l

reader = pd.read_csv(lab_path, iterator=True)
lab_df = pd.DataFrame(columns=lab_cols)

# N > 839995358
n = 0

tmp = reader.get_chunk(30000000)
tmp = tmp[tmp.OMOP_ID.isin([3050449,3049506,3020876,3020564])]
tmp = tmp[lab_cols]
lab_df = pd.concat([lab_df,tmp])
lab_df.to_csv('df_lab_values.csv', index=None)

# lab_df = pd.read_csv('df_lab_values.csv')
# ep_df = pd.read_csv('df_nephropathy.csv')

##############################################################
# Prepare well-defined endpoints and demographic information
##############################################################

first_event_path = '/data/processed_data/endpointer/archive/wide_first_events_DF10_2022_09_29.txt.ALL.gz'
df_events = pd.read_csv(first_event_path, sep='\t', usecols=['FINREGISTRYID','T1D','T2D','FG_CVD'])
# ,'DM1_NEPHROPATHY'
df_events.to_csv('df_dm_cvd.csv', index=None)

info_path = '/data/processed_data/minimal_phenotype/archive/minimal_phenotype_2022-03-28.csv'
df_info = pd.read_csv(info_path)
df = df_info[['FINREGISTRYID','id_father','id_mother','sex','date_of_birth','death_date','emigration_date']]
df = df[df.FINREGISTRYID.isin(df_events[df_events.T1D == 1].FINREGISTRYID)]

################################################
# Combine all the information collected above
################################################

# nothing like E102&N083, only one code in one cell
# len(ep_df[ep_df['E102&N083'] == True]) = 0

# this includes all those diagnosed as nephropathy at least once
# N = 1618 among 36439
df['nephropathy_ep1'] = np.select([(~df.FINREGISTRYID.isin(ep_df.FINREGISTRYID)),
                               (df.FINREGISTRYID.isin(ep_df.FINREGISTRYID))],[0.0,1.0])

df = df[~df.sex.isna()] # total: 36439 -> 36152; cases: 1618 -> 1581

# select those with at least two diagnoses
tmp = ep_df.FINREGISTRYID.value_counts()
tmp = tmp[(tmp > 1) == True]
ep_df = ep_df[ep_df.FINREGISTRYID.isin(tmp.keys())] # 58071

# N = 1069 among 36439
df['nephropathy_ep2'] = np.select([(~df.FINREGISTRYID.isin(ep_df.FINREGISTRYID)),
                               (df.FINREGISTRYID.isin(ep_df.FINREGISTRYID))],[0.0,1.0])

ep_df = ep_df.sort_values(['FINREGISTRYID','EVENT_DAY'])
ep_df = ep_df[ep_df.duplicated(subset=['FINREGISTRYID']) == False]
df = df.merge(ep_df[['FINREGISTRYID','EVENT_AGE','EVENT_DAY']], 'left')
df = df.rename(columns={'EVENT_AGE':'age_ep','EVENT_DAY':'day_ep'})

