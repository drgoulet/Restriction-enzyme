# Search for restriction enzyme cut sites within a given amino acid sequence
# Input: Amino acid query sequence
# Output: Restriction enzymes and where they cut in the translated DNA sequence
# Dennis R. Goulet
# First upload to Github: 17 May 2022

import itertools

aa_seq = "GARC"

enzymes = {
	"AatII":"GACGTC",
	"AclI":"AACGTT",
	"AfeI":"AGCGCT",
	"AflII":"CTTAAG",
	"AgeI":"ACCGGT",
	"ApaI":"GGGCCC",
	"ApaLI":"GTGCAC",
	"AscI":"GGCGCGCC",
	"AseI":"ATTAAT",
	"AsiSI":"GCGATCGC",
	"AvrII":"CCTAGG",
	"BamHI":"GGATCC",
	"BclI":"TGATCA",
	"BglII":"AGATCT",
	"BsiWI":"CGTACG",
	"BspEI":"TCCGGA",
	"BspHI":"TCATGA",
	"BsrGI":"TGTACA",
	"BssHII":"GCGCGC",
	"BstBI":"TTCGAA",
	"BstZ17I":"GTATAC",
	"ClaI":"ATCGAT",
	"DraI":"TTTAAA",
	"EagI":"CGGCCG",
	"EcoRI":"GAATTC",
	"EcoRV":"GATATC",
	"FseI":"GGCCGGCC",
	"FspI":"TGCGCA",
	"HindIII":"AAGCTT",
	"HpaI":"GTTAAC",
	"KasI":"GGCGCC",
	"KpnI":"GGTACC",
	"MfeI":"CAATTG",
	"MluI":"ACGCGT",
	"MscI":"TGGCCA",
	"NcoI":"CCATGG",
	"NdeI":"CATATG",
	"NheI":"GCTAGC",
	"NotI":"GCGGCCGC",
	"NruI":"TCGCGA",
	"NsiI":"ATGCAT",
	"PacI":"TTAATTAA",
	"PciI":"ACATGT",
	"PmeI":"GTTTAAAC",
	"PmlI":"CACGTG",
	"PsiI":"TTATAA",
	"PstI":"CTGCAG",
	"PvuI":"CGATCG",
	"PvuII":"CAGCTG",
	"SacI":"GAGCTC",
	"SalI":"GTCGAC",
	"SbfI":"CCTGCAGG",
	"ScaI":"AGTACT",
	"SnaBI":"TACGTA",
	"SpeI":"ACTAGT",
	"SphI":"GCATGC",
	"SrfI":"GCCCGGGC",
	"SspI":"AATATT",
	"StuI":"AGGCCT",
	"SwaI":"ATTTAAAT",
	"XbaI":"TCTAGA",
	"XhoI":"CTCGAG",
	"XmaI":"CCCGGG"
}

codon_table = {
'A': ('GCT', 'GCC', 'GCA', 'GCG'),
'C': ('TGT', 'TGC'),
'D': ('GAT', 'GAC'),
'E': ('GAA', 'GAG'),
'F': ('TTT', 'TTC'),
'G': ('GGT', 'GGC', 'GGA', 'GGG'),
'H': ('CAT', 'CAC'),
'I': ('ATT', 'ATC', 'ATA'),
'K': ('AAA', 'AAG'),
'L': ('TTA', 'TTG', 'CTT', 'CTC', 'CTA', 'CTG'),
'M': ('ATG',),
'N': ('AAT', 'AAC'),
'P': ('CCT', 'CCC', 'CCA', 'CCG'),
'Q': ('CAA', 'CAG'),
'R': ('CGT', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'),
'S': ('TCT', 'TCC', 'TCA', 'TCG', 'AGT', 'AGC'),
'T': ('ACT', 'ACC', 'ACA', 'ACG'),
'V': ('GTT', 'GTC', 'GTA', 'GTG'),
'W': ('TGG',),
'Y': ('TAT', 'TAC'),
}

dna_seq = []
index = 0

for aa in aa_seq:
	for codon in codon_table[aa]:
		dna_seq.append(codon)

dna_list = list(itertools.product(*[list(codon_table[key]) for key in aa_seq]))
dna_list_flat = [''.join(a) for a in dna_list ]

blank_seq = []
for j in range(0,3*len(aa_seq)):
	blank_seq.append('-')
tuple(blank_seq)
blank_str = ''.join(blank_seq)
final_list = []
new_seq_str = ''

for seq in dna_list_flat:
	for enzyme, recog in enzymes.items():
		result = seq.find(recog)
		if result == -1:
			pass
		else:
			new_seq = blank_seq[:]
			for i in range(0,len(seq)):
					if seq[i:i+len(recog)] == recog:
						res = i
						break
			for k in range(res,res+len(recog)):
				new_seq[k] = seq[k]
			new_seq_str = ''.join(new_seq)

			if new_seq_str in final_list:
				pass
			else:
				print(f"The enzyme {enzyme} cuts over '{recog}' in the sequence {new_seq_str}.")
				final_list.append(new_seq_str)
				new_seq_str = ''
