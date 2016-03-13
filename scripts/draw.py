from Bio import Phylo

tree = Phylo.read("./phyloxml.xml", "phyloxml")
print(tree)
Phylo.draw(tree)
