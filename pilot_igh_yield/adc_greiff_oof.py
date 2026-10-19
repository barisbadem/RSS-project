"""Out-of-frame yield of the Greiff 2017 mouse IGH repertoires (ADC study PRJEB18631) on VDJServer.

usage: python3 -I adc_greiff_oof.py mouse_repertoires.json out.tsv

Same rule as the earlier human survey: frame comes from junction_length (vj_in_frame is null on
VDJServer), and a junction must be >= 15 nt to count (productive=false is often a failed annotation).
"""
import json
import subprocess
import sys
import time

HOST = "https://vdjserver.org"
MIN_JUNCTION = 15


def post(endpoint, body):
    proc = subprocess.run(
        ["curl", "-sS", "--max-time", "400", "-H", "Content-Type: application/json",
         "-d", json.dumps(body), f"{HOST}/airr/v1/{endpoint}"],
        capture_output=True, text=True)
    try:
        return json.loads(proc.stdout)
    except Exception:
        return None


def facet(rep_id, field, extra=None):
    content = [{"op": "=", "content": {"field": "repertoire_id", "value": rep_id}}]
    if extra:
        content.append(extra)
    flt = content[0] if len(content) == 1 else {"op": "and", "content": content}
    body = {"filters": flt, "facets": field}
    for attempt in range(6):
        out = post("rearrangement", body)
        if out and "Facet" in out:
            return out["Facet"]
        time.sleep(6 * (attempt + 1))
    return None


def main():
    data = json.load(open(sys.argv[1]))
    reps = [r for r in data["vdjserver"] if r["study"]["study_id"] == "PRJEB18631"
            and ((r["sample"][0].get("cell_subset") or {}).get("label") in ("precursor B cell", "naive B cell"))]
    rows = []
    for r in reps:
        s, smp = r["subject"], r["sample"][0]
        cohort = smp.get("disease_state_sample")
        subset = (smp.get("cell_subset") or {}).get("label")
        prod = facet(r["repertoire_id"], "productive")
        npj = facet(r["repertoire_id"], "junction_length",
                    {"op": "=", "content": {"field": "productive", "value": False}})
        if prod is None or npj is None:
            print("FAILED", s["subject_id"], flush=True)
            continue
        n_prod = sum(f["count"] for f in prod if f.get("productive") is True)
        n_np = sum(f["count"] for f in prod if f.get("productive") is False)
        oof = sum(f["count"] for f in npj if f.get("junction_length") and f["junction_length"] >= MIN_JUNCTION and f["junction_length"] % 3)
        rows.append((subset, cohort, s["subject_id"], s.get("strain_name"), n_prod, n_np, oof))
        print(subset, cohort, s["subject_id"], n_prod, n_np, oof, flush=True)
    with open(sys.argv[2], "w") as fh:
        fh.write("cell_subset\tcohort\tsubject\tstrain\tproductive\tnon_productive_raw\tout_of_frame_ge15nt\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
