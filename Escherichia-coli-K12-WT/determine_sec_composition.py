from Bio import SeqIO
from collections import Counter

with open('secretion_machinery.fasta', 'rU') as fin:
	for record in SeqIO.parse(fin, 'fasta'):
		print record.id
		cnt = Counter(str(record.seq))
		print '''
		      <composition>
        <componentReference component="A" stoichiometry="{}"/>
        <componentReference component="C" stoichiometry="{}"/>
        <componentReference component="E" stoichiometry="{}"/>
        <componentReference component="D" stoichiometry="{}"/>
        <componentReference component="G" stoichiometry="{}"/>
        <componentReference component="F" stoichiometry="{}"/>
        <componentReference component="I" stoichiometry="{}"/>
        <componentReference component="H" stoichiometry="{}"/>
        <componentReference component="K" stoichiometry="{}"/>
        <componentReference component="M" stoichiometry="{}"/>
        <componentReference component="L" stoichiometry="{}"/>
        <componentReference component="N" stoichiometry="{}"/>
        <componentReference component="Q" stoichiometry="{}"/>
        <componentReference component="P" stoichiometry="{}"/>
        <componentReference component="S" stoichiometry="{}"/>
        <componentReference component="R" stoichiometry="{}"/>
        <componentReference component="T" stoichiometry="{}"/>
        <componentReference component="W" stoichiometry="{}"/>
        <componentReference component="V" stoichiometry="{}"/>
        <componentReference component="Y" stoichiometry="{}"/>
      </composition>
      '''.format(cnt['A'],
      			 cnt['C'],
      			 cnt['E'],
      			 cnt['D'],
      			 cnt['G'],
      			 cnt['F'],
      			 cnt['I'],
      			 cnt['H'],
      			 cnt['K'],
      			 cnt['M'],
      			 cnt['L'],
      			 cnt['N'],
      			 cnt['Q'],
      			 cnt['P'],
      			 cnt['S'],
      			 cnt['R'],
      			 cnt['T'],
      			 cnt['W'],
      			 cnt['V'],
      			 cnt['Y'],)

