"""Turn two truncated .fastq.gz mate files into plain FASTQ with the same number of
complete, correctly paired records.

usage: python3 -I slice_pair.py R1.partial.gz R2.partial.gz out_R1.fq out_R2.fq [max_pairs]

A byte-range download cuts the gzip stream mid-record, and the two mates compress
differently, so they end at different reads. This keeps only records that are
complete in both files and whose read ids (first token of the header) agree.
Prints the number of pairs written.
"""
import gzip
import sys
import zlib


def records(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        try:
            while True:
                head = fh.readline()
                if not head:
                    return
                seq = fh.readline()
                plus = fh.readline()
                qual = fh.readline()
                if not qual.endswith("\n") or len(qual) != len(seq):
                    return  # truncated final record
                yield head, seq, plus, qual
        except (EOFError, zlib.error, OSError):
            return


def main() -> None:
    r1_path, r2_path, out1, out2 = sys.argv[1:5]
    max_pairs = int(sys.argv[5]) if len(sys.argv) > 5 else None
    n = 0
    with open(out1, "w") as o1, open(out2, "w") as o2:
        for a, b in zip(records(r1_path), records(r2_path)):
            if a[0].split()[0] != b[0].split()[0]:
                print(f"read id mismatch at pair {n + 1}: {a[0].strip()} vs {b[0].strip()}", file=sys.stderr)
                break
            o1.write("".join(a))
            o2.write("".join(b))
            n += 1
            if max_pairs and n >= max_pairs:
                break
    print(n)


if __name__ == "__main__":
    main()
