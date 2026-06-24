Scripts for telomeric variant repeat (TVR) and motif analysis of long-read nanopore sequencing data, comparing ALT+ and ALT− samples.
All scripts read/write paths under `/omics/groups/OE0436/data/linmq/...` on the DKFZ cluster — edit the hardcoded paths before reuse.

## Cluster job submission

`submit_script.sh`
LSF batch job: loads Python/minimap2, installs dependencies, calls `run_script.sh`.

`run_script.sh`
Runs Telogator2 on one SRR sample with ONT settings, producing `tlens_by_allele.tsv` and `tel_reads.fa.gz`.

## TVR / motif counting

`tvr_count.py`
Counts TVR symbols in `tvr_consensus` at offsets from the sequence start; writes `ALT+_all_counts.tsv`.

`reverse tvr_count.py`
Same as above but trims from the end; near-duplicate of `tvr_count.py` aside from trim direction.

`matched.py`
Trims reads to their read_TL cutoff and scans for `keys.json` motifs (fwd + reverse complement); writes normalized motif counts per SRR.

`matched to extract raw reads.py`
Companion to `matched.py`; extracts the matching read sequences themselves instead of counting motifs.

`Extracting Chr from tsv`
Filters reads by chromosome + allele and extracts them from the gzipped FASTA (no file extension, but it's Python — called `extract.py` in usage notes).

## Sequence manipulation / simulation

`Simulation- taking out Consensus Sequences`
Filters by chromosome + allele and writes `tvr_consensus` sequences out as FASTA.

`reduction.py`
Greedy 6-mer matcher that reduces sequences to single-letter TVR codes, choosing the better-reducing orientation (forward vs reverse complement).

`combination for simulated sequences.py`
Builds synthetic test sequences from fixed telomere repeat fragments flanked by `TTAGGG` repeats.

`percentage.py`
Meant to append a percentage of one sequence onto another. Has a bug: the line defining `portion` is commented out, so it currently raises a `NameError`.

## Plotting

`DFHIN.py`
Boxplots of TVRs D/F/H/I/N by offset, ALT+ vs ALT−.

`plot according to tvr.py`
Three-group (ALT+/ALT−/negative control) line plots per TVR with SEM shading.

`top tvr.py`
Top-5 TVR bar charts at a fixed offset, ALT+ vs ALT−.

`tvr-stacked.py`
Clustered stacked bar chart of TVR composition by cell type, plus a stats table.

`plot-matched.py` / `plot by motifs.py`
Plot `matched.py`'s motif tables: the first as a forward/backward bar chart, the second as a batch stripplot over multiple motif lengths.

## Reference data & notes

`kmers.tsv`
K-mer lookup table used by the motif-counting scripts.

`Reading_files.txt`
Scratchpad of example CLI commands showing how the scripts chain together end to end (some referenced names don't match files in the repo exactly).
