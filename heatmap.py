import os
import sys
import argparse
import subprocess
from scipy import stats

parser = argparse.ArgumentParser('Create a matrix that can be loaded into R and generate a heatmap of a metric')

parser.add_argument('gtex_compare_dir', help='Directory containing all comparisons made. Folders are named tissue1-vs-tissue2')
parser.add_argument('ccle_compare_dir', help='Directory containing all comparisons made. Folders are named tissue1-vs-tissue2')
#parser.add_argument('-m', '--metrics', nargs='+', default=('skew'), choices=('mean', 'median', 'skew', 'mode', 'ratio', 'mannwhitneyu', 'delta'), help='A metric to generate the heatmap for')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to')

args = parser.parse_args()

if not os.path.isdir(args.outdir):
    try:
        os.makedirs(args.outdir)
    except OSError:
        pass

def dicToMtx(dic, equal=1):
    rownames = colnames = dic.keys()
    length = len(rownames)
    header = '\t' + ('\t').join(rownames) + '\n'
    out = [header]
    for i in xrange(length):
        ki = rownames[i]
        line = ['{}'.format(ki)]
        for j in xrange(length):
            kj = rownames[j]
            if ki == kj:
                line.append(equal)
                continue
            try:
                line.append(dic[ki][kj])
            except KeyError:
                line.append(dic[kj][ki])
        line = ('\t').join([str(x) for x in line])
        out.append(line)
    return ('\n').join(out)

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

distance_mtx = {}

for compare in os.listdir(args.ccle_compare_dir):
    distances_ccle = os.path.join(args.ccle_compare_dir, compare, 'distances')
    distances_gtex = os.path.join(args.gtex_compare_dir, compare, 'distances')
    distances_ccle = parseDistances(distances_ccle)
    distances_gtex = parseDistances(distances_gtex)
    pval = stats.mannwhitneyu(distances_ccle, distances_gtex)[1]
    t1, t2 = compare.split('-vs-')
    if t1 not in distance_mtx:
        distance_mtx[t1] = {t2: pval}
    elif t2 not in distance_mtx[t1]:
        distance_mtx[t1][t2] = pval

result = dicToMtx(distance_mtx, equal=0)

out = os.path.join(args.outdir, 'delta.mtx')
with open(out, 'w') as o:
    o.write(result)

rscript = '/gsc/btl/linuxbrew/bin/Rscript'
script = os.path.join(args.outdir, 'script.R')
with open(script, 'w') as o:
    o.write('library(pheatmap)\n')
    o.write('library(corrplot)\n')
    o.write('setwd(\'{}\')\n'.format(args.outdir))
    o.write('data_delta = read.table(\'delta.mtx\')\n')
    o.write('png(\'correlation_delta.png\', width=1000, height=1100, res=100)\n')
    o.write('corrplot(cor(data_delta), method="number", type="upper", mar = c(4,4,4,4))\n')
    o.write('dev.off()\n')
    o.write('png(\'heatmap_delta.png\', width=1000, height=750, res=150)\n')
    o.write('pheatmap(data_delta, cluster_rows=FALSE, cluster_cols=FALSE)\n')
    o.write('dev.off()\n')

subprocess.check_call([rscript, script])
