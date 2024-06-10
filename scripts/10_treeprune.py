from ete3 import Tree
import re
import os
treefile = open("/gpfs/scratch/withomas/eve/data/trees/4705sp_mean.nwk","r")
tree = treefile.read()
t = Tree(tree, format = 0)
t.prune(["sorex_araneus","rattus_norvegicus","mus_musculus","cavia_aperea","bos_taurus","sus_scrofa","homo_sapiens","macaca_mulatta","pan_troglodytes","fukomys_mechowii","capra_hircus","ovis_aries","ictidomys_tridecemlineatus","peromyscus_maniculatus","microtus_ochrogaster","heterocephalus_glaber"])
print(t)
t.write(format = 1, outfile = "hypot_ALVAREZ_TreePruned.nh")

#dropout
treefile = open("/gpfs/scratch/withomas/eve/data/trees/4705sp_mean.nwk","r")
tree = treefile.read()
t = Tree(tree, format = 0)
t.prune(["rattus_norvegicus","mus_musculus","cavia_aperea","bos_taurus","sus_scrofa","homo_sapiens","macaca_mulatta","pan_troglodytes","fukomys_mechowii","capra_hircus","ovis_aries","ictidomys_tridecemlineatus","peromyscus_maniculatus","microtus_ochrogaster","heterocephalus_glaber"])
print(t)
t.write(format = 1, outfile = "hypot_ALVAREZ_TreePruned_DO.nh")
