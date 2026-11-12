#!/usr/bin/env python3
"""
Architecture of construct reads, and the integration loci they define.

Reads `samtools view` output on stdin, WITH supplementary alignments:

  samtools view -F 0x104 subset/<id>.construct_reads.bam \
    | architecture.py --construct <contig> [--construct-fa <construct.fa>]

Each read is one independent observation of one molecule, so nothing aligns
reads to each other -- the structure is already in the alignments.

WHAT IT REPORTS
  per read     the ordered segments, as [construct_start-end strand] and
               <genome_contig:pos strand>
  loci         construct<->genome junctions grouped by position. A real
               integration site is FIXED on the junction side and SCATTERED on
               the shear side; distinct shear points count independent
               molecules, which duplicates and duplex read-pairs cannot fake.
  array        construct bases sitting between the junction and the read end --
               how many copies the molecule actually carries above the junction
  junctions    construct-construct joins: orientation, breakpoints, and whether
               they fall in regions where alignment cannot be trusted

UNRELIABLE REGIONS (--construct-fa auto-detects, or give --mask)
  Terminal inverted repeats: both construct ends carry the same sequence, so a
  segment lying inside one cannot be placed to an end, and single-stranded DNA
  snaps back there during library prep. Both inflate apparent inverted joins.
  Long homopolymers: ONT undercounts them, so alignments carry a deletion or
  terminate there, producing breakpoints that are basecalling, not biology.
  Breakpoints in these regions are FLAGGED, never dropped -- for a linear
  replicon the telomere is exactly where real integration junctions sit.
"""
import sys, re, argparse
from collections import Counter, defaultdict

CIG = re.compile(r'(\d+)([MIDNSHP=X])')
RC = str.maketrans('ACGTacgt', 'TGCAtgca')


def parse_cigar(c):
    ops = [(int(n), o) for n, o in CIG.findall(c)]
    lead = ops[0][0] if ops and ops[0][1] in 'SH' else 0
    tail = ops[-1][0] if len(ops) > 1 and ops[-1][1] in 'SH' else 0
    qaln = sum(n for n, o in ops if o in 'MI=X')
    rspan = sum(n for n, o in ops if o in 'MDN=X')
    return lead, qaln, tail, rspan


def find_unreliable(path):
    """Terminal inverted repeats and long homopolymers in the construct."""
    seq = ''.join(l.strip() for l in open(path) if not l.startswith('>')).upper()
    n = len(seq)
    out = []

    # terminal inverted repeat: longest prefix that is the revcomp of the suffix
    rc = seq.translate(RC)[::-1]
    lo, hi = 0, n // 2
    while lo < hi:                                   # binary search on length
        mid = (lo + hi + 1) // 2
        if seq[:mid] == rc[:mid]:
            lo = mid
        else:
            hi = mid - 1
    if lo >= 20:
        out.append((1, lo, f'TIR {lo}bp'))
        out.append((n - lo + 1, n, f'TIR {lo}bp'))

    # homopolymer runs
    i = 0
    while i < n:
        j = i
        while j + 1 < n and seq[j + 1] == seq[i]:
            j += 1
        if j - i + 1 >= 12:
            out.append((i + 1, j + 1, f'{j-i+1}x{seq[i]}'))
        i = j + 1
    return n, out


def flag(pos, unreliable):
    for a, b, label in unreliable:
        if a <= pos <= b:
            return label
    return None


def cluster(values, within):
    """-> list of (modal value, members) for values grouped within `within`."""
    out, cur = [], []
    for v in sorted(values):
        if cur and v - cur[-1] > within:
            out.append(cur); cur = []
        cur.append(v)
    if cur:
        out.append(cur)
    return out


def load_segments(stream, min_bp=100, only=None):
    """Parse `samtools view` output into {read name: [segment, ...]}.

    entry/exit are the reference coordinates where the read enters and leaves a
    segment, which is what junction logic needs: on a minus-strand alignment the
    read runs from the HIGH coordinate to the low one.
    """
    segs = defaultdict(list)
    for line in stream:
        f = line.rstrip('\n').split('\t')
        if len(f) < 6 or f[0].startswith('@'):
            continue
        name, flg, ref, pos, mapq, cig = f[0], int(f[1]), f[2], int(f[3]), int(f[4]), f[5]
        if only and name not in only:
            continue
        if flg & 0x4 or flg & 0x100:          # unmapped or secondary
            continue
        lead, qaln, tail, rspan = parse_cigar(cig)
        if qaln < min_bp:
            continue
        rev = bool(flg & 0x10)
        rstart = (tail if rev else lead) + 1   # minus-strand clips are on the revcomp read
        rs, re_ = pos, pos + rspan - 1
        segs[name].append(dict(
            read_start=rstart, read_end=rstart + qaln - 1, qaln=qaln,
            ref=ref, ref_start=rs, ref_end=re_, strand='-' if rev else '+',
            mapq=mapq, total=lead + qaln + tail,
            entry=(re_ if rev else rs),
            exit=(rs if rev else re_)))
    for v in segs.values():
        v.sort(key=lambda d: d['read_start'])
    return segs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--construct', required=True)
    ap.add_argument('--construct-fa', help='auto-detect TIRs and homopolymers')
    ap.add_argument('--mask', help='extra unreliable regions, e.g. 1-110,1985-2035')
    ap.add_argument('--min-bp', type=int, default=100)
    ap.add_argument('--reads', help='restrict to these read names (comma separated)')
    ap.add_argument('--locus-within', type=int, default=500,
                    help='junctions this close are one locus (default 500)')
    ap.add_argument('--max-reads-shown', type=int, default=40)
    args = ap.parse_args()

    clen, unreliable = (0, [])
    if args.construct_fa:
        clen, unreliable = find_unreliable(args.construct_fa)
    if args.mask:
        for part in args.mask.split(','):
            a, b = part.split('-')
            unreliable.append((int(a), int(b), 'masked'))

    only = set(args.reads.split(',')) if args.reads else None
    C = args.construct
    segs = load_segments(sys.stdin, min_bp=args.min_bp, only=only)

    if not segs:
        sys.exit('no alignments read on stdin')

    if unreliable:
        print(f'construct {C}' + (f' ({clen} bp)' if clen else ''))
        print('unreliable regions (breakpoints here are flagged *):')
        for a, b, label in sorted(set(unreliable)):
            print(f'  {a:>6}-{b:<6} {label}')
        print()

    # Genome junction positions, binned, so each can be called recurrent or a
    # one-off. A chimera formed during adapter ligation joins two sequences that
    # each map perfectly well -- MAPQ says nothing about whether the JOIN is
    # real. Only recurrence across independently sheared molecules does.
    BIN = 50_000
    recur = Counter()
    for s_ in segs.values():
        for a, b in zip(s_, s_[1:]):
            if a['ref'] == C and b['ref'] != C:
                recur[(b['ref'], b['entry'] // BIN)] += 1
            elif a['ref'] != C and b['ref'] == C:
                recur[(a['ref'], a['exit'] // BIN)] += 1

    def is_recurrent(ref, pos):
        return recur[(ref, pos // BIN)] >= 3

    # ---------------- per read ----------------
    loci = defaultdict(list)
    copies_hist, orient, arch_counter = Counter(), Counter(), Counter()
    cc_breaks = Counter()

    order = sorted(segs.items(), key=lambda kv: -kv[1][0]['total'])
    print(f'{"read":<40} {"len":>7}  architecture')
    print('-' * 112)
    for i, (name, s) in enumerate(order):
        ncop = sum(1 for d in s if d['ref'] == C)
        copies_hist[ncop] += 1
        arch_counter[''.join('C' if d['ref'] == C else 'G' for d in s)] += 1

        if i < args.max_reads_shown:
            parts = []
            for d in s:
                if d['ref'] == C:
                    parts.append(f"[{d['ref_start']}-{d['ref_end']}{d['strand']}]")
                else:
                    mark = '**' if is_recurrent(d['ref'], d['entry']) else ''
                    parts.append(f"<{d['ref']}:{d['ref_start']}{d['strand']}>{mark}")
            print(f'{name[:40]:<40} {s[0]["total"]:>7}  ' + ' '.join(parts))

        # construct bases carried by this molecule, and how many unit-lengths
        cbp = sum(d['qaln'] for d in s if d['ref'] == C)

        for a, b in zip(s, s[1:]):
            if a['ref'] == C and b['ref'] == C:
                orient['direct' if a['strand'] == b['strand'] else 'INVERTED'] += 1
                cc_breaks[(a['exit'], b['entry'])] += 1
            elif a['ref'] == C and b['ref'] != C:
                loci[b['ref']].append(dict(read=name, gj=b['entry'], shear=b['exit'],
                                           cj=a['exit'], cbp=cbp, strand=b['strand'],
                                           glen=b['qaln']))
            elif a['ref'] != C and b['ref'] == C:
                loci[a['ref']].append(dict(read=name, gj=a['exit'], shear=a['entry'],
                                           cj=b['entry'], cbp=cbp, strand=a['strand'],
                                           glen=a['qaln']))

    if len(order) > args.max_reads_shown:
        print(f'... and {len(order)-args.max_reads_shown} more reads')

    # ---------------- integration loci ----------------
    print()
    allj = [(contig, j) for contig, v in loci.items() for j in v]
    rec_j = [j for c, j in allj if is_recurrent(c, j['gj'])]
    one = [j for c, j in allj if not is_recurrent(c, j['gj'])]

    def med(xs, key):
        v = sorted(k[key] for k in xs)
        return v[len(v) // 2] if v else 0

    print()
    print('=' * 112)
    print('JUNCTION RECURRENCE   (** marks a recurrent genome anchor above)')
    print(f'  recurrent {len(rec_j):>4}   median genome-side anchor '
          f'{med(rec_j, "glen"):,} bp')
    print(f'  one-off   {len(one):>4}   median genome-side anchor '
          f'{med(one, "glen"):,} bp'
          f'   ({sum(1 for j in one if j["glen"] < 500)} of them under 500 bp)')
    print()
    print('  A one-off sits at a random construct AND a random genome coordinate, so it')
    print('  is not evidence of integration whatever its MAPQ -- MAPQ says the piece is')
    print('  placed unambiguously, not that the JOIN is real. The short ones are mostly')
    print('  spurious brief alignments; the long ones look like ligation chimeras. Only')
    print('  recurrence across independently sheared molecules separates real from not.')
    print()
    print('=' * 112)
    print('INTEGRATION LOCI  (junction fixed + shear scattered = real; '
          'distinct shear points = independent molecules)')
    print()
    rows = []
    for contig, js in loci.items():
        for grp in cluster([j['gj'] for j in js], args.locus_within):
            lo, hi = min(grp), max(grp)
            mem = [j for j in js if lo <= j['gj'] <= hi]
            shears = {j['shear'] for j in mem}
            cjs = Counter(j['cj'] for j in mem)
            rows.append((len(mem), contig, lo, hi, shears, cjs, mem))
    rows.sort(reverse=True, key=lambda r: (r[0], len(r[4])))

    if not rows:
        print('  none -- no construct<->genome junctions')
    for nread, contig, lo, hi, shears, cjs, mem in rows[:15]:
        span = max(shears) - min(shears)
        cj_mode, cj_n = cjs.most_common(1)[0]
        fl = flag(cj_mode, unreliable)
        print(f'  {contig}:{lo}' + (f'-{hi}' if hi != lo else ''))
        print(f'      reads {nread}   independent molecules {len(shears)}   '
              f'junction spread {hi-lo} bp   shear spread {span:,} bp')
        print(f'      construct side joins at base {cj_mode}'
              f' ({cj_n}/{nread} reads agree)' + (f'  * {fl}' if fl else ''))
        arr = sorted(j['cbp'] for j in mem)
        if clen:
            units = [f'{a/clen:.1f}' for a in arr]
            print(f'      construct carried per molecule (unit lengths): '
                  f'{", ".join(units[:12])}' + (' ...' if len(units) > 12 else ''))
        if len(shears) < nread:
            print(f'      NOTE only {len(shears)} distinct shear points for {nread} reads'
                  f' -- possible duplicates or duplex pairs')
        print()

    # ---------------- copy number ----------------
    print('=' * 112)
    print('CONSTRUCT COPIES PER READ')
    for k in sorted(copies_hist):
        print(f'  {k:>2} segments  {copies_hist[k]:>4} reads  {"#" * min(copies_hist[k], 70)}')

    # ---------------- construct-construct junctions ----------------
    tot = sum(orient.values())
    if tot:
        print()
        print('=' * 112)
        print('CONSTRUCT-CONSTRUCT JUNCTIONS')
        for k, v in orient.most_common():
            print(f'  {k:<10} {v:>5}  ({100*v/tot:.0f}%)')
        print()
        print('  most common breakpoint pairs (exit -> entry):')
        for (x, e), n in cc_breaks.most_common(12):
            fx, fe = flag(x, unreliable), flag(e, unreliable)
            mark = ' *' if (fx or fe) else ''
            why = f'   [{fx or "-"} -> {fe or "-"}]' if mark else ''
            print(f'    {x:>6} -> {e:<6}  {n:>4}{mark}{why}')
        flagged = sum(n for (x, e), n in cc_breaks.items()
                      if flag(x, unreliable) or flag(e, unreliable))
        print(f'\n  {flagged}/{tot} ({100*flagged/tot:.0f}%) of these junctions fall in an '
              f'unreliable region -- treat copy number and orientation with caution')

    print()
    print('architecture signatures (C=construct, G=genome):')
    for a, n in arch_counter.most_common(10):
        print(f'  {a:<26} {n:>4}')


if __name__ == '__main__':
    main()
