import os
import sys
import argparse
import subprocess

parser = argparse.ArgumentParser('Create a matrix that can be loaded into R and generate a heatmap of a metric')

parser.add_argument('gtex_compare_dir', help='Directory containing all comparisons made. Folders are named tissue1-vs-tissue2')
parser.add_argument('ccle_compare_dir', help='Directory containing all comparisons made. Folders are named tissue1-vs-tissue2')
parser.add_argument('-m', '--metrics', nargs='+', default=('skew'), choices=('mean', 'median', 'skew', 'mode', 'ratio', 'mannwhitneyu', 'delta'), help='A metric to generate the heatmap for')
#parser.add_argument('-c', '--create', action='store_true', help='Whether or not to run R and generate a heatmap')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to')

args = parser.parse_args()

if not os.path.isdir(args.outdir):
    try:
        os.makedirs(args.outdir)
    except OSError:
        pass

def parseDistances(_file):
    distances = []
    with open(_file, 'r') as f:
        for line in f:
            line = float(line.strip())
            distances.append(line)
    return distances

def parseStats(path):
    res = {}
    with open(path, 'r') as f:
        for line in f:
            line = line.strip().split('\t')
            line[0] = line[0][:-1]
            res[line[0]] = float(line[1])
    return res

results_gtex = {}

for compare in os.listdir(args.gtex_compare_dir):
    stats = os.path.join(args.gtex_compare_dir, compare, 'statistics')
    #distances = os.path.join(args.gtex_compare_dir, compare, 'distances')
    mwu = os.path.join(args.gtex_compare_dir, compare, 'mannwhitneyu')
    stats = parseStats(stats)
    #distances = parseDistances(distances)
    stats['distances'] = distances
    stats['mannwhitneyu'] = None
    if 'mannwhitneyu' in args.metrics:
        stats['mannwhitneyu'] = parseStats(mwu)['mannwhitneyu']
    stats['ratio'] = abs(stats['total_pos'] / stats['total_neg'])
    t1, t2 = compare.split('-vs-')
    if t1 not in results_gtex:
        results_gtex[t1] = {t2: stats}
    elif t2 not in results_gtex[t1]:
        results_gtex[t1][t2] = stats

results_ccle = {}

for compare in os.listdir(args.ccle_compare_dir):
    stats = os.path.join(args.ccle_compare_dir, compare, 'statistics')
    #distances = os.path.join(args.gtex_compare_dir, compare, 'distances')
    stats = parseStats(stats)
    mwu = os.path.join(args.ccle_compare_dir, compare, 'mannwhitneyu')
    stats['mannwhitneyu'] = None
    #stats['distances'] = distances
    if 'mannwhitneyu' in args.metrics:
        stats['mannwhitneyu'] = parseStats(mwu)['mannwhitneyu']
    stats['ratio'] = abs(float(stats['total_pos']) / stats['total_neg'])
    t1, t2 = compare.split('-vs-')
    if t1 not in results_ccle:
        results_ccle[t1] = {t2: stats}
    elif t2 not in results_ccle[t1]:
        results_ccle[t1][t2] = stats

keys = sorted(results_ccle.keys())
len_keys = len(keys)
for metric in args.metrics:
    out = os.path.join(args.outdir, metric + '.mtx')
    with open(out, 'w') as o:
        o.write('\t' + ('\t').join(keys) + '\n')
        for i,t1 in enumerate(keys):
            results = results_ccle
            o.write(t1)
            for j,t2 in enumerate(keys):
                if t1 == t2:
                    if metric == 'ratio':
                        o.write('\t1')
                    else:
                        o.write('\t0')
                    results = results_gtex
                    continue
                try:
                    o.write('\t')
                    o.write(str(results[t1][t2][metric]))
                except KeyError:
                    try:
                        o.write('\t')
                        o.write(str(results[t2][t1][metric]))
                    except KeyError:
                        o.write('\t0')
            o.write('\n')

rscript = '/gsc/btl/linuxbrew/bin/Rscript'
script = os.path.join(args.outdir, 'script.R')
with open(script, 'w') as o:
    o.write('library(pheatmap)\n')
    o.write('library(corrplot)\n')
    o.write('setwd(\'{}\')\n'.format(args.outdir))
    for metric in args.metrics:
        o.write('data_{} = read.table(\'{}.mtx\')\n'.format(metric, metric))
        o.write('png(\'correlation_{}.png\', width=1000, height=1100, res=100)\n'.format(metric))
        o.write('corrplot(cor(data_{}), method="number", type="upper", mar = c(4,4,4,4))\n'.format(metric))
        o.write('dev.off()\n')
        o.write('png(\'heatmap_{}.png\', width=1000, height=750, res=150)\n'.format(metric))
        o.write('pheatmap(data_{}, cluster_rows=FALSE, cluster_cols=FALSE)\n'.format(metric))
        o.write('dev.off()\n')

subprocess.check_call([rscript, script])
