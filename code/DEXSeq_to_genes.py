#! /usr/bin/env python2

'''
Converts exon expression of DEXSeq results to average gene expression values
for a GO enrichment analysis.

# input.csv:


# output.csv:


#command:

$  python DEXSeq_to_GOinput.py \
  -i input.csv \
  -c log2fold \
  -o output.csv

#contact:

Dmytro Kryvokhyzha dmytro.kryvokhyzha@evobio.eu

'''

############################# modules #############################

import argparse, sys  # for input options

############################# options #############################

class CommandLineParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

parser = CommandLineParser()
parser.add_argument(
    '-i',
    '--input',
    help='name of the DEXSeq results table',
    type=str,
    required=True)
parser.add_argument(
    '-c',
    '--colname',
    help='name of the DEXSeq results column to average',
    type=str,
    required=True)
parser.add_argument(
    '-o',
    '--output',
    help='name of the output file',
    type=str, required=True)
args = parser.parse_args()


############################# program #############################

colname = args.colname
output = open(args.output, 'w')
geneDict = {}

with open(args.input) as datafile:

  header = datafile.readline().strip('\n').replace('"', '').split("\t")
  try:
    colindex = header.index(colname)
  except ValueError:
    print('Invalid column name %s' % colname)

  output.write('gene\tmean_%s\n' % str(header[colindex]))

  for line in datafile:
    words = line.strip('\n').replace('"', '').split("\t")
    group = words[0]
    if words[colindex] != "NA":
      logFC = float(words[colindex])
      for i in group.split('+'):
        # split overlapping exons:
        gene = i.split('_')[0]
        # create a dict of genes and their logFC:
        if gene in geneDict.keys(): 
          geneDict[gene].append(logFC)
        else:
          geneDict[gene] = [logFC]

for k in geneDict.keys():
  l = geneDict[k]
  avg_logFC = sum(l) / len(l)
  output.write('%s\t%s\n' % (k, avg_logFC))

datafile.close()
output.close()

print('Done!')
