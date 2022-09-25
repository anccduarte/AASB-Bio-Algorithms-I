
# -*- coding: utf-8 -*-

import re

class Enzymes:
	"""Verificação dos locais de corte e dos fragmentos originados pela ação de uma enzima de restrição."""

	def __init__(self, seq: str, enzyme: str) -> None:

		if type(seq) != str or type(enzyme) != str:
			raise TypeError("Os parâmetros 'seq' e 'enzyme' devem ser do tipo 'string'.")

		for c in seq.upper():
			if c not in "ACTG":
				raise ValueError("A sequência que inseriu não corresponde a DNA.")

		for c in enzyme.upper():
			if c not in "^ACTGRYMKSWBDHVN":
				raise ValueError("Por favor, insira uma enzima de restrição válida.")

		self.seq = seq.upper()
		self.enzyme = enzyme.upper()


	def __str__(self) -> str:
		"""Imprime os fragmentos de DNA originados pelo corte da enzima de uma forma legível."""
		return ("Fragmentos de DNA:\n" + " | ".join(self.cut_subseqs()))

	
	def __iub2dna(self, iub: str) -> str:
		"""Recebe uma sequência de caracteres em formato IUPAC-IUB e devolve uma sequência de caracteres correspondentes a nucleótidos."""
		poss = "R [GA] Y [CT] M [AC] K [GT] S [GC] W [AT] B [^A] D [^C] H [^G] V [^T] N [ATCG]".split()
		dic = {poss[i]: poss[i+1] for i in range(0,len(poss),2)}
		iub = list(iub)
		for i,c in enumerate(iub):
			if c in dic.keys():
				iub[i] = dic[c]
		return "".join(iub)


	def __enzyme2re(self) -> str:
		"""Recebe uma sequência de caracteres pertencentes a uma enzima e devolve a expressão regular correspondente."""
		one,two = self.enzyme.split("^")
		one = self.__iub2dna(one)
		two = self.__iub2dna(two)
		return f"(?<={one}){two}"


	def cut_positions(self) -> list:
		"""Recebe uma enzima e uma sequência de DNA e devolve uma lista de índices correspondentes às posições de corte."""
		enz = self.__enzyme2re()
		pattern = re.compile(enz)
		matches = pattern.finditer(self.seq)
		return [match.start() for match in matches]


	def cut_subseqs(self) -> list:
		"""Recebe uma enzima e uma sequência de DNA e devolve uma lista de fragmentos originados pelo corte da enzima."""
		cuts = [0]
		for cut in self.cut_positions():
			cuts.append(cut)
		return [self.seq[cuts[i]:cuts[i+1]] if i+1 < len(cuts) else self.seq[cuts[i]:] for i in range(len(cuts))]

