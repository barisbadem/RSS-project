"""Inventory of mouse (Mus musculus) IGH repertoires exposed by AIRR Data Commons repositories.

usage: python3 -I adc_mouse_inventory.py out.json
Writes the raw repertoire metadata; summarising is done in a second step.
"""
import json
import subprocess
import sys

HOSTS = {
    "vdjserver": "https://vdjserver.org",
    "ipa1": "https://ipa1.ireceptor.org",
    "ipa2": "https://ipa2.ireceptor.org",
    "ipa3": "https://ipa3.ireceptor.org",
    "ipa4": "https://ipa4.ireceptor.org",
    "ipa5": "https://ipa5.ireceptor.org",
    "ipa6": "https://ipa6.ireceptor.org",
}


def post(host, endpoint, body):
    proc = subprocess.run(
        ["curl", "-sS", "--max-time", "300", "-H", "Content-Type: application/json",
         "-d", json.dumps(body), f"{HOSTS[host]}/airr/v1/{endpoint}"],
        capture_output=True, text=True,
    )
    try:
        return json.loads(proc.stdout)
    except Exception:
        return {"_error": proc.stdout[:300] + proc.stderr[:300]}


def main():
    out = {}
    for host in HOSTS:
        reps, start = [], 0
        while True:
            body = {
                "filters": {"op": "=", "content": {"field": "subject.species.label", "value": "Mus musculus"}},
                "size": 1000,
                "from": start,
            }
            res = post(host, "repertoire", body)
            page = res.get("Repertoire")
            if page is None:
                print(host, "ERROR", str(res)[:200])
                break
            reps.extend(page)
            if len(page) < 1000:
                break
            start += 1000
        out[host] = reps
        print(host, len(reps), "mouse repertoires")
    with open(sys.argv[1], "w") as fh:
        json.dump(out, fh)


if __name__ == "__main__":
    main()
