
# -*- coding: utf-8 -*-

class Seq:
	"""Implementação de métodos básicos para a análise de sequências de DNA."""

	def __init__(self, seq: str) -> None:

		if type(seq) != str:
			raise TypeError("A sequência deve ser do tipo 'string'.")

		for c in seq:
			if c.upper() not in "ACTG":
				raise ValueError("A sequência que inseriu não corresponde a DNA.")

		self.seq = seq.upper()


	def __str__(self) -> str:
		"""Imprime o código utilizado aquando da criação da instância."""
		return f'Seq(seq = "{self.seq}")'


	def get_codons(self) -> list:
		"""Recebe uma sequência de DNA e devolve uma lista dos respetivos codões."""
		rna = self.seq.replace("T","U")
		return [rna[i:i+3] for i in range(0,len(rna),3) if i+3 <= len(rna)]


	def get_amino(self) -> str:
		"""Recebe uma sequência de DNA e devolve uma sequência de aminoácidos."""
		codons = self.get_codons()
		dic = {"F": ["UUU","UUC"], "L": ["UUA","UUG","CUU","CUC","CUA","CUG"], "I": ["AUU","AUC","AUA"], "M": ["AUG"],
			   "V": ["GUU","GUC","GUA","GUG"], "S": ["UCU","UCC","UCA","UCG","AGU","AGC"], "P": ["CCU","CCC","CCA","CCG"],
			   "T": ["ACU","ACC","ACA","ACG"], "A": ["GCU","GCC","GCA","GCG"], "Y": ["UAU","UAC"], "_": ["UAA","UAG","UGA"],
			   "H": ["CAU","CAC"], "Q": ["CAA","CAG"], "N": ["AAU","AAC"], "K": ["AAA","AAG"], "D": ["GAU","GAC"], "E": ["GAA","GAG"],
			   "C": ["UGU","UGC"], "W": ["UGG"], "R": ["CGU","CGC","CGA","CGG","AGA","AGG"], "G": ["GGU","GGC","GGA","GGG"]}
		amino = ""
		for codon in codons:
			for k in dic:
				if codon in dic[k]:
					amino += k
		return amino


	def get_orfs(self) -> list:
		"""Recebe uma sequência de DNA e devolve uma lista com as seis ORFS."""
		dna = self.seq
		inv = dna[::-1]
		comp_inv = inv.lower().replace("a","T").replace("t","A").replace("c","G").replace("g","C")
		return [dna,dna[1:],dna[2:],comp_inv,comp_inv[1:],comp_inv[2:]]


	def get_prots(self) -> list:
		"""Recebe uma sequência de DNA e devolve uma lista das possíveis proteínas presentes na reading frame corrente."""
		prots_list = []
		aminos = self.get_amino()
		for i,amino in enumerate(aminos):
			if amino == "M":
				prot = ""
				for j in range(i,len(aminos)):
					if aminos[j] != "_":
						prot += aminos[j]
					else:
						break
				prots_list.append(prot)
		return prots_list


	def get_all_prots(self) -> list:
		"""Recebe uma sequência de DNA e devolve uma lista com todas as proteínas possíveis."""
		orfs = self.get_orfs()
		all_prots = []
		for orf in orfs:
			orf = Seq(orf)
			prots = orf.get_prots()
			all_prots.append(prots)
		return [prots for nested_prots in all_prots for prots in nested_prots]

