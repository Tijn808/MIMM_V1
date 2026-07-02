#!/usr/bin/env python3
"""
cohort_table1.py  --  Build Table 1 as the 15-region validation set of Sisman et
al.: the 10 JHU white-matter ROIs plus the 5 subcortical grey-matter structures,
MIMM MVF vs T2-GRASE MWF (mean MVF, mean MWF, bias = MVF - MWF, Cohen's d).

It just merges the two CSVs already produced upstream:
  cohort_mwf.py     -> cohort_mwf_paper_wm.csv   (10 WM ROIs, L+R combined)
  cohort_gm_rois.py -> cohort_gm_rois.csv        (5 GM structures)

Run those two first, then this. Output:
  <results>/cohort_analysis/cohort_table1_15region.csv

Usage: python3 cohort_table1.py <results_dir>
"""
import sys, os, csv

if len(sys.argv) < 2:
    sys.exit('usage: cohort_table1.py <results_dir>')
ca = os.path.join(sys.argv[1], 'cohort_analysis')
wm_csv = os.path.join(ca, 'cohort_mwf_paper_wm.csv')
gm_csv = os.path.join(ca, 'cohort_gm_rois.csv')
for p in (wm_csv, gm_csv):
    if not os.path.exists(p):
        sys.exit(f'missing {p}\n'
                 'Run cohort_mwf.py and cohort_gm_rois.py first, then rerun this.')


def rows_from(path, name_col):
    out = []
    with open(path) as f:
        for r in csv.DictReader(f):
            out.append([r[name_col], r['MVF_mean'], r['MWF_mean'],
                        r['bias_MVF_minus_MWF'], r['cohens_d']])
    return out


wm = rows_from(wm_csv, 'tract')
gm = rows_from(gm_csv, 'structure')

out = os.path.join(ca, 'cohort_table1_15region.csv')
with open(out, 'w', newline='') as f:
    w = csv.writer(f)
    w.writerow(['Region', 'MVF', 'MWF', 'Bias (MVF-MWF)', 'Cohen d'])
    w.writerow(['White matter (JHU)', '', '', '', ''])
    w.writerows(wm)
    w.writerow(['Subcortical grey matter (Harvard-Oxford)', '', '', '', ''])
    w.writerows(gm)

print(f'saved: {out}   ({len(wm)} WM + {len(gm)} GM = {len(wm) + len(gm)} regions)')
print(f'{"Region":46s} {"MVF":>7s} {"MWF":>7s} {"Bias":>8s} {"d":>6s}')
for r in wm + gm:
    print(f'{r[0]:46s} {r[1]:>7s} {r[2]:>7s} {r[3]:>8s} {r[4]:>6s}')
