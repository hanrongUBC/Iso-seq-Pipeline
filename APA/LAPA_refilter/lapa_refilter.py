#!/usr/bin/env python3
"""Re-filter a finished LAPA run with explicit per-cluster rules.

LAPA does not need to be re-run. Starting from the pooled clusters saved
before any sample, internal-priming or replication filtering
(results/raw_polyA_clusters.bed), this script:

  1. recovers integer per-sample counts from LAPA's same-strand raw count
     BigWigs (results/counts/<sample>_polyA_counts_{pos,neg}.bw);
  2. scores internal priming from the genome sequence downstream of the
     cleavage site, case-insensitively and without LAPA's exemption for
     clusters that carry a human poly(A)-signal hexamer;
  3. applies an explicit replicate rule inside each treatment group,
     replacing LAPA's 1,000-cluster rolling replication rate;
  4. writes one count matrix (union of kept clusters x all samples) for a
     replicate-aware test such as run_dexseq_apa.R.

Coordinates: LAPA polyA_site is a 0-based boundary. The last transcribed
base is polyA_site (1-based) on + and polyA_site + 1 (1-based) on -.
LAPA's own BED Start/End put minus-strand ends one base outside the read;
the site_bed_start/site_bed_end columns written here are corrected.
"""
import argparse
import csv
import gzip
import json
from collections import OrderedDict
from pathlib import Path

import pyBigWig

COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')


class Fasta:
    """Minimal indexed FASTA reader (needs <fasta>.fai); returns uppercase."""

    def __init__(self, path):
        self.handle = open(path, 'rb')
        self.index = {}
        with open(str(path) + '.fai') as fai:
            for line in fai:
                name, length, offset, line_bases, line_bytes = line.split('\t')[:5]
                self.index[name] = tuple(map(int, (length, offset, line_bases, line_bytes)))

    def fetch(self, chrom, start, end):
        length, offset, line_bases, line_bytes = self.index[chrom]
        start, end = max(0, start), min(length, end)
        if start >= end:
            return ''
        first = offset + start // line_bases * line_bytes + start % line_bases
        last = offset + (end - 1) // line_bases * line_bytes + (end - 1) % line_bases + 1
        self.handle.seek(first)
        raw = self.handle.read(last - first).decode('ascii')
        return raw.replace('\n', '').replace('\r', '').upper()


def downstream(fasta, chrom, strand, site, size):
    """Genomic bases just after the cleavage site, in transcript orientation."""
    if strand == '+':
        return fasta.fetch(chrom, site, site + size)
    return fasta.fetch(chrom, site - size, site).translate(COMPLEMENT)[::-1]


def longest_a_run(seq):
    best = run = 0
    for base in seq:
        run = run + 1 if base == 'A' else 0
        best = max(best, run)
    return best


def read_samples(path):
    with open(path) as handle:
        rows = list(csv.DictReader(handle, delimiter='\t'))
    if not rows or not {'sample', 'group'} <= set(rows[0]):
        raise ValueError('samples table needs tab-separated columns: sample, group')
    groups = OrderedDict()
    for row in rows:
        groups.setdefault(row['group'], []).append(row['sample'])
    samples = [row['sample'] for row in rows]
    if len(set(samples)) != len(samples):
        raise ValueError('duplicate sample names in samples table')
    return samples, groups


def read_clusters(path):
    """Group rows by (chrom, start, end, strand); a cluster overlapping
    several genes appears once with all its gene_ids."""
    clusters = OrderedDict()
    with open(path) as handle:
        for n, row in enumerate(csv.reader(handle, delimiter='\t'), 1):
            if len(row) < 14:
                raise ValueError(f'{path}:{n}: expected 14 LAPA columns, found {len(row)}')
            chrom, start, end, site, _count, strand, feature, gene = row[:8]
            if strand not in ('+', '-'):
                raise ValueError(f'{path}:{n}: unstranded cluster')
            key = (chrom, int(start), int(end), strand)
            entry = clusters.setdefault(key, {'site': int(site), 'features': [], 'genes': []})
            if entry['site'] != int(site):
                raise ValueError(f'{path}:{n}: conflicting polyA_site for {key}')
            if gene not in entry['genes']:
                entry['genes'].append(gene)
                entry['features'].append(feature)
    return clusters


def count_interval(bw, chroms, chrom, start, end):
    if chrom not in chroms:
        return 0
    spans = bw.intervals(chrom, start, end) or ()
    total = sum((min(b, end) - max(a, start)) * v for a, b, v in spans)
    if total < 0 or abs(total - round(total)) > 1e-5:
        raise ValueError('expected raw integer read-end counts, not normalised BigWigs')
    return int(round(total))


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--lapa-results', required=True, type=Path, help='LAPA results/ directory')
    p.add_argument('--fasta', required=True, type=Path,
                   help='the same genome FASTA LAPA used (indexed, .fai next to it)')
    p.add_argument('--samples', required=True, type=Path, help='TSV with columns sample, group')
    p.add_argument('--out', required=True, type=Path)
    p.add_argument('--control', default='Control')
    p.add_argument('--min-count', type=int, default=3,
                   help='reads a sample needs at a cluster to count as observed (LAPA used 1)')
    p.add_argument('--min-reps', type=int, default=2,
                   help='observed replicates needed inside at least one group')
    p.add_argument('--a10', type=int, default=7, help='internal priming if >= this many A in 10 nt downstream')
    p.add_argument('--a20', type=int, default=12, help='internal priming if >= this many A in 20 nt downstream')
    p.add_argument('--a-run', type=int, default=6, help='internal priming if an A run this long occurs in 20 nt')
    args = p.parse_args()

    samples, groups = read_samples(args.samples)
    if args.control not in groups:
        raise ValueError(f'control group {args.control!r} not in samples table')
    clusters = read_clusters(args.lapa_results / 'raw_polyA_clusters.bed')
    final_keys = set(read_clusters(args.lapa_results / 'polyA_clusters.bed'))
    keys = list(clusters)
    fasta = Fasta(args.fasta)
    missing = sorted({k[0] for k in keys} - set(fasta.index))
    if missing:
        raise ValueError(f'cluster contigs absent from FASTA index (wrong FASTA?): {missing[:5]}')

    counts = {key: [0] * len(samples) for key in keys}
    for j, sample in enumerate(samples):
        for strand, suffix in (('+', 'pos'), ('-', 'neg')):
            path = args.lapa_results / 'counts' / f'{sample}_polyA_counts_{suffix}.bw'
            with pyBigWig.open(str(path)) as bw:
                chroms = bw.chroms()
                for key in keys:
                    if key[3] == strand:
                        counts[key][j] = count_interval(bw, chroms, key[0], key[1], key[2])
        print(f'counts recovered: {sample}', flush=True)

    index = {s: i for i, s in enumerate(samples)}
    args.out.mkdir(parents=True, exist_ok=True)
    stats = OrderedDict((name, 0) for name in (
        'raw_clusters', 'internal_priming', 'replicated', 'kept', 'kept_annotated_gene',
        'lapa_final', 'kept_and_lapa_final', 'kept_not_in_lapa_final', 'lapa_final_not_kept',
        'multi_gene_clusters'))
    header = (['cluster_id', 'chrom', 'start', 'end', 'strand', 'polyA_site', 'last_base_1based',
               'site_bed_start', 'site_bed_end', 'feature', 'gene_id', 'n_genes',
               'A_in_10nt', 'A_in_20nt', 'longest_A_run', 'internal_priming']
              + [f'{s}' for s in samples]
              + [f'observed_reps_{g}' for g in groups] + ['replicated', 'in_lapa_final', 'keep'])
    kept_rows = []
    with gzip.open(args.out / 'clusters_refiltered.tsv.gz', 'wt') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        for i, key in enumerate(keys, 1):
            chrom, start, end, strand = key
            entry = clusters[key]
            site = entry['site']
            seq20 = downstream(fasta, chrom, strand, site, 20)
            a10, a20, run = seq20[:10].count('A'), seq20.count('A'), longest_a_run(seq20)
            primed = a10 >= args.a10 or a20 >= args.a20 or run >= args.a_run
            row_counts = counts[key]
            reps = {g: sum(row_counts[index[s]] >= args.min_count for s in members)
                    for g, members in groups.items()}
            replicated = any(n >= args.min_reps for n in reps.values())
            keep = replicated and not primed
            in_final = key in final_keys
            genes = [g for g in entry['genes'] if g and g != '-1' and not g.startswith('intergenic_')]
            gene = ';'.join(genes) if genes else entry['genes'][0]
            feature = ';'.join(entry['features'])
            cid = f'C{i:07d}'
            # Last transcribed base as a 0-based BED interval: [site-1, site) on +, [site, site+1) on -.
            bed_start = site - 1 if strand == '+' else site
            last_base = site if strand == '+' else site + 1
            writer.writerow([cid, chrom, start, end, strand, site, last_base, bed_start, bed_start + 1,
                             feature, gene, len(genes), a10, a20, run, int(primed)]
                            + row_counts + [reps[g] for g in groups]
                            + [int(replicated), int(in_final), int(keep)])
            stats['raw_clusters'] += 1
            stats['internal_priming'] += primed
            stats['replicated'] += replicated
            stats['kept'] += keep
            stats['multi_gene_clusters'] += len(genes) > 1
            stats['lapa_final'] += in_final
            stats['kept_and_lapa_final'] += keep and in_final
            stats['kept_not_in_lapa_final'] += keep and not in_final
            stats['lapa_final_not_kept'] += in_final and not keep
            if keep:
                stats['kept_annotated_gene'] += len(genes) == 1
                kept_rows.append([cid, chrom, strand, site, last_base, feature, gene, len(genes)] + row_counts)

    with open(args.out / 'count_matrix.tsv', 'w') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        writer.writerow(['cluster_id', 'chrom', 'strand', 'polyA_site', 'last_base_1based',
                         'feature', 'gene_id', 'n_genes'] + samples)
        writer.writerows(kept_rows)

    summary = {'counts': stats,
               'parameters': {'min_count': args.min_count, 'min_reps': args.min_reps,
                              'internal_priming': {'A_in_10nt_at_least': args.a10,
                                                   'A_in_20nt_at_least': args.a20,
                                                   'A_run_at_least': args.a_run,
                                                   'signal_hexamer_exemption': False,
                                                   'case_insensitive': True}},
               'groups': groups, 'source': str(args.lapa_results.resolve())}
    (args.out / 'refilter_summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main()
