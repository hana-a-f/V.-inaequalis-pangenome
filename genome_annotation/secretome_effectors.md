# Effector Prediction

This folder documents the effector-candidate prediction pipeline. Each isolate's representative protein set is filtered down to a candidate secretome — proteins with a signal peptide **and** no transmembrane domain — and that secretome is then scanned with **EffectorP 3.0** to flag candidate effectors.

## Pipeline overview

1. **SignalP 5.0** — predict signal peptides (identifies proteins routed for secretion)
2. Custom script — pull SignalP-positive proteins into a FASTA
3. **TMHMM 2.0c** — predict transmembrane helices
4. Custom script — filter TMHMM output to those with ≥1 predicted TM helix
5. Custom script — extract those TM-containing proteins into a FASTA
6. Custom script — subtract TM proteins from SignalP proteins → **secretome**
7. Collapse to one representative per locus
8. **EffectorP 3.0** — score each secreted protein as effector candidate / non-effector
9. Strip EffectorP boilerplate from the output for downstream parsing

## Logic

A **secreted protein** = has a signal peptide AND no transmembrane domain. Signal peptide alone isn't enough — a protein anchored in the membrane also gets a signal peptide but never actually leaves the cell. Subtracting TMHMM-positive proteins from SignalP-positive ones removes those false positives and leaves genuine candidate secreted proteins for EffectorP to score.

## Tools used

| Tool | Version / source |
| --- | --- |
| SignalP | 5.0 (`/programs/signalp-5.0/bin/signalp`) |
| TMHMM | 2.0c (`/programs/TMHMM2.0c/bin/`) |
| EffectorP | 3.0 (from https://github.com/JanaSperschneider/EffectorP-3.0) |
| Weka | 3-8-4 (bundled with EffectorP, unzipped after clone) |
| GNU parallel | (system) |

## Custom scripts you need to add to this folder

**⚠ Pull these four Python scripts from your working directory on the cluster and upload them into this folder:**

| Script | Purpose |
| --- | --- |
| `get_signalp.py` | Takes a SignalP `.gff3` + the query protein FASTA, writes a FASTA of only the proteins with a predicted signal peptide |
| `filter_predicted_hel.py` | Takes a raw TMHMM `.tmhmm.out`, writes a filtered version keeping only entries with predicted transmembrane helices |
| `make_transmembrane_protein_fasta.py` | Takes the filtered TMHMM output + the query protein FASTA, writes a FASTA of only the TM-containing proteins |
| `secretome.py` | Takes SignalP-positive FASTA + TM-positive FASTA + output path, writes a FASTA of proteins that pass SignalP but don't have TM domains — the candidate secretome |



## Step 1 — SignalP 5.0

```bash
/programs/signalp-5.0/bin/signalp --version   # sanity check
mkdir -p signalp
> signalp.cmds

for f in *.fa; do
    iso="${f%.fa}"
    echo "/programs/signalp-5.0/bin/signalp \
  -org euk \
  -batch 20000 \
  -fasta $f \
  -format short \
  -gff3 \
  -prefix signalp/${iso}" >> signalp.cmds
done

parallel -j 2 < signalp.cmds
```

Key flags: `-org euk` (eukaryotic), `-batch 20000` (batch size for performance), `-gff3` (also write a GFF3 of hits — used in Step 2).

## Step 2 — Pull SignalP-positive proteins into a FASTA

```bash
touch signalp.py.cmds

for f in signalp/*.gff3; do
    iso="${f%%.gff3}"
    prot="${iso}.fa"
    echo "python get_signalp.py $f $prot ${iso}.signalp.fa" >> signalp.py.cmds
done

parallel -j 2 < signalp.py.cmds
```

## Step 3 — TMHMM 2.0c

```bash
export PATH=/programs/TMHMM2.0c/bin:$PATH
touch tmhmm.cmds

for f in *.fa; do
    iso="${f%%.fa}"
    echo "tmhmm --short < \"$f\" > ${iso}.tmhmm.out" >> tmhmm.cmds
done

parallel -j 5 < tmhmm.cmds
```

## Step 4 — Filter TMHMM output

Keeps only proteins with predicted transmembrane helices:

```bash
touch filter.tmhmm.cmds

for f in *.tmhmm.out; do
    iso="${f%%.tmhmm.out}"
    out="${iso}.tmhmm.filtered.out"
    echo "python filter_predicted_hel.py $f $out" >> filter.tmhmm.cmds
done

parallel -j 5 < filter.tmhmm.cmds
```

## Step 5 — Pull TM-containing proteins into a FASTA

```bash
touch transmembrane.cmds

for f in *.fa; do
    iso="${f%%.fa}"
    in="${iso}.tmhmm.filtered.out"
    out="${iso}.transmembrane.proteins.fa"
    echo "python make_transmembrane_protein_fasta.py $in $f $out" >> transmembrane.cmds
done

parallel -j 5 < transmembrane.cmds
```

## Step 6 — Secretome = SignalP MINUS TM

```bash
touch secretome.cmds

for f in *.signalp.fa; do
    iso="${f%.signalp.fa}"
    t="${iso}.transmembrane.proteins.fa"
    out="${iso}.secretome.fa"
    echo "python secretome.py $f $t $out" >> secretome.cmds
done

parallel -j 5 < secretome.cmds
```

## Step 7 — Filter secretome to one gene per locus

Keep only one isoform per locus (strips `.t<N>` transcript suffix) — matches the collapsed-by-locus convention from `genome_annotation/`:

```bash
#!/bin/bash
set -euo pipefail

for PEPTIDE in VI_*.fa; do
    ISO=$(basename "$PEPTIDE" .fa)
    SECRETOME="${ISO}.secretome.fa"
    OUT="${ISO}.locus.filt.secretome.fa"

    [[ -s "$SECRETOME" ]] || continue
    echo "Processing $ISO"

    # Extract locus IDs from secretome
    grep "^>" "$SECRETOME" \
      | sed 's/^>//' \
      | cut -d' ' -f1 \
      | sed 's/\.t[0-9]\+$//' \
      | sort -u \
      > "${ISO}.secretome.loci.txt"

    # Filter peptide FASTA
    awk '
    BEGIN {
        while ((getline < "'"${ISO}.secretome.loci.txt"'") > 0)
            keep[$1]=1
    }
    /^>/ {
        id=$1
        sub(/^>/,"",id)
        sub(/\.t[0-9]+$/,"",id)
        print_flag = (id in keep)
    }
    print_flag { print }
    ' "$PEPTIDE" > "$OUT"

    rm "${ISO}.secretome.loci.txt"
done

echo "=== ALL ISOLATES DONE ==="
```

## Step 8 — Install and run EffectorP 3.0

```bash
git clone https://github.com/JanaSperschneider/EffectorP-3.0.git
cd EffectorP-3.0/
unzip weka-3-8-4.zip
cd ..

touch effectp.cmds

for f in *.locus.filt.secretome.fa; do
    iso="${f%%.locus.filt.secretome.fa}"
    echo "python EffectorP.py -i $f > ${iso}.effectp.txt" >> effectp.cmds
done

parallel -j 10 < effectp.cmds
```

## Step 9 — Clean up EffectorP output for parsing

Strip the header banner (first 10 lines) and footer (last 13 lines) to leave only the tabular result rows:

```bash
for f in *.effectp.txt; do
    base="${f%%.effectp.txt}"
    out="${base}.effectp.filtered.txt"
    tail -n +11 "$f" | head -n -13 > "$out"
done
```

**Note:** the `+11` / `-13` line counts are specific to the exact EffectorP 3.0 output format seen when we ran it. If you rerun with a different EffectorP version and the banner/footer size changes, these numbers will silently cut real data or leave boilerplate in. A safer version would use a text marker (e.g. `awk '/^#/{f=0} /Identifier/{f=1} f'`), but the numbers work for the current run.

## Outputs

Per isolate:
- `<iso>.secretome.fa` — full candidate secretome (before locus filtering)
- `<iso>.locus.filt.secretome.fa` — locus-collapsed secretome (input to EffectorP)
- `<iso>.effectp.txt` — raw EffectorP output
- `<iso>.effectp.filtered.txt` — cleaned tabular EffectorP output for parsing
