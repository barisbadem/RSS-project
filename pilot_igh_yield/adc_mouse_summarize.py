import json, sys
from collections import defaultdict, Counter

data = json.load(open(sys.argv[1]))

def first(x):
    if isinstance(x, list):
        return x[0] if x else {}
    return x or {}

def label(d):
    d = first(d)
    return (d.get("label") if isinstance(d, dict) else None) or ""

studies = defaultdict(list)
for host, reps in data.items():
    for r in reps:
        st = r.get("study", {})
        studies[(host, st.get("study_id"), st.get("study_title"))].append(r)

print(f"{len(studies)} studies\n")
for (host, sid, title), reps in sorted(studies.items(), key=lambda kv: -len(kv[1])):
    strains, dx, tissue, subset, locus, tmpl, subj = Counter(), Counter(), Counter(), Counter(), Counter(), Counter(), set()
    for r in reps:
        s = r.get("subject", {}); subj.add(s.get("subject_id"))
        strains[s.get("strain_name") or "?"] += 1
        d = first(s.get("diagnosis")); dx[label(d.get("disease_diagnosis")) or "none/healthy?"] += 1
        for smp in (r.get("sample") or []):
            tissue[label(smp.get("tissue")) or "?"] += 1
            subset[label(smp.get("cell_subset")) or "?"] += 1
            tmpl[smp.get("template_class") or "?"] += 1
            for pt in (smp.get("pcr_target") or []):
                locus[pt.get("pcr_target_locus") or "?"] += 1
    print(f"[{host}] {sid} | n_rep={len(reps)} n_subj={len(subj)}")
    print(f"   title : {(title or '')[:140]}")
    print(f"   strain: {dict(strains.most_common(4))}")
    print(f"   dx    : {dict(dx.most_common(5))}")
    print(f"   tissue: {dict(tissue.most_common(4))} | subset: {dict(subset.most_common(4))}")
    print(f"   locus : {dict(locus.most_common(4))} | template: {dict(tmpl.most_common(3))}\n")
