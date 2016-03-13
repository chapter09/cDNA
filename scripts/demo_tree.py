from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor
from StringIO import StringIO
from Bio.Align.Applications import MuscleCommandline
from Bio import AlignIO
from Bio import Phylo
from Bio import SeqIO
from Bio import Seq
import os

from Bio.Align.Applications import ClustalwCommandline


#input_file = '../data/GCF_000008805.1_ASM880v1_protein.faa'
input_file = '../data/GCF_000008805.1_ASM880v1_protein.faa'
#records = SeqIO.parse(input_file, 'fasta')
#records = list(records) # make a copy, otherwise our generator
#                        # is exhausted after calculating maxlen
#maxlen = max(len(record.seq) for record in records)

## pad sequences so that they all have the same length
#for record in records:
#    if len(record.seq) != maxlen:
#        sequence = str(record.seq).ljust(maxlen, '.')
#        record.seq = Seq.Seq(sequence)
#assert all(len(record.seq) == maxlen for record in records)

## write to temporary file and do alignment
#output_file = '{}_padded.fasta'.format(os.path.splitext(input_file)[0])
#with open(output_file, 'w') as f:
#    SeqIO.write(records, f, 'fasta')
#alignment = AlignIO.read(output_file, "fasta")

#cline = ClustalwCommandline("clustalw2", infile=input_file)
#print(cline)
#print type(cline)

muscle_cline = MuscleCommandline(input=input_file)
stdout, stderr = muscle_cline()
alignment = AlignIO.read(StringIO(stdout), "fasta")
print(alignment)

#alignment = AlignIO.read('../data/ls_orchid.fasta', 'fasta')
#print alignment
calculator = DistanceCalculator('ident')
dm = calculator.get_distance(alignment)
constructor = DistanceTreeConstructor()
tree = constructor.upgma(dm)
Phylo.write(tree, 'phyloxml.xml', 'phyloxml')
