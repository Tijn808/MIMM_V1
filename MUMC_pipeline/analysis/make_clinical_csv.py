#!/usr/bin/env python3
"""
make_clinical_csv.py  —  turn the Castor IMPROMYMS export into a clean clinical.csv
for cohort_clinical.py.

The export is semicolon-delimited with a BOM. The only field that needs decoding is
the EDSS: t1_edss_score is stored as an option CODE, and real EDSS = code * 0.5
(code 0 -> 0.0, 2 -> 1.0, 4 -> 2.0, 11 -> 5.5, ...). Everything else (SDMT count,
mean T25ftWT, 9-Hole-Peg trial times) is already numeric.

Outputs columns:
  participant_id, age, sex, bmi, edss, sdmt, t25ftwt, hpt_dom, hpt_nondom
Rows with no EDSS date (healthy controls / not-yet-assessed) keep NaN in the
clinical columns; the correlation script masks those out.

Usage:  python3 make_clinical_csv.py <IMPROMYMS_export.csv> [clinical.csv]
        (default output: clinical.csv next to this script)

NOTE: clinical.csv contains pseudonymised patient data — do not commit it to git.
"""
import sys, os
import pandas as pd

if len(sys.argv) < 2:
    sys.exit('usage: make_clinical_csv.py <IMPROMYMS_export.csv> [clinical.csv]')
src = sys.argv[1]
out = sys.argv[2] if len(sys.argv) > 2 else 'clinical.csv'

df = pd.read_csv(src, sep=';', encoding='utf-8-sig', dtype=str)
df.columns = [c.strip() for c in df.columns]


def num(col):
    """numeric column or all-NaN if the column is absent."""
    if col not in df.columns:
        return pd.Series([float('nan')] * len(df))
    return pd.to_numeric(df[col].str.replace(',', '.', regex=False), errors='coerce')


hpt_dom = pd.concat([num('t1_dominant_trial_1_time'),
                     num('t1_dominant_hand_trial_2_time')], axis=1).mean(axis=1).round(2)
hpt_nondom = pd.concat([num('t1_non_dominant_trial_1_time'),
                        num('t1_non_dominant_trial_2_time')], axis=1).mean(axis=1).round(2)

clin = pd.DataFrame({
    'participant_id': df['participant_id'].str.strip(),
    'age':       num('age_baseline'),
    'sex':       num('sex'),                 # 0 = man, 1 = woman (Castor coding)
    'bmi':       num('t0_bmi'),
    'edss':      num('t1_edss_score') * 0.5,  # decode option code -> real EDSS
    'sdmt':      num('t1_number_sdmt'),
    't25ftwt':   num('t1_mean_t25ftwt'),
    'hpt_dom':   hpt_dom,
    'hpt_nondom': hpt_nondom,
})

clin.to_csv(out, index=False)
n_edss = clin['edss'].notna().sum()
n_sdmt = clin['sdmt'].notna().sum()
print(f'wrote {out}: {len(clin)} participants '
      f'({n_edss} with EDSS, {n_sdmt} with SDMT)')
if n_edss:
    e = clin['edss'].dropna()
    print(f'EDSS  range {e.min():.1f}-{e.max():.1f}, median {e.median():.1f}')
if n_sdmt:
    s = clin['sdmt'].dropna()
    print(f'SDMT  range {s.min():.0f}-{s.max():.0f}, median {s.median():.0f}')
