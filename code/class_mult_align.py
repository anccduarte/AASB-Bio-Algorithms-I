
# -*- coding: utf-8 -*-

class MultipleAlign:
	"""Construção de alinhamentos múltiplos (alinhamento progressivo)."""

	def __init__(self, seqs: list, scoring = [2,0], gap = 4) -> None:

		if type(seqs) != list:
			raise TypeError("As sequências a alinhar devem estar contidas numa lista.")

		for seq in seqs:
			if type(seq) != str:
				raise TypeError("As sequências que compõem o parâmetro 'seqs' devem ser do tipo 'str'.")

		for seq in seqs:
			for c in seq:
				if c.upper() not in "ACTG":
					raise ValueError("Pelo menos uma das sequências não corresponde a DNA.")

		if type(scoring) != list:
			raise TypeError("O parâmetro 'scoring' deve ser do tipo 'list'.")

		if type(scoring[0]) != int or type(scoring[1]) != int:
			raise TypeError("Os valores que compõem o parâmetro 'scoring' devem ser do tipo 'int'.")

		if type(gap) != int:
			raise TypeError("O parâmetro 'gap' deve ser do tipo 'int'.")

		if gap < 0:
			raise ValueError("O valor do parâmetro 'gap' não deve ser menor que 0.")

		self.seqs = []
		for seq in seqs:
			self.seqs.append(seq.upper())
		self.scoring = scoring
		self.gap = gap


	def __str__(self) -> str:
		"""Devolve a lista de sequências original e o alinhamento múltiplo impresso de uma forma mais legível."""
		align = "\n".join(self.mult_align())
		return f"Sequences: {self.seqs}\nMultiple Alignment:\n{align}"


	def __calc_score(self, x1: str, x2: str) -> int:
		"""Recebe dois caracteres e devolve o score correspondente."""
		return self.scoring[0] if x1 == x2 else self.scoring[1]


	def __build_mats(self, seq1: str, seq2: str) -> list:
		"""Recebe 2 sequências e o valor da gap penalty, e retorna as matrizes de score e trace."""

		nrows = len(seq1) + 1
		ncols = len(seq2) + 1
		score_mat = [[0 for j in range(ncols)] for i in range(nrows)]
		trace_mat = [[0 for j in range(ncols)] for i in range(nrows)]

		for i in range(1,nrows):
			score_mat[i][0] = score_mat[i-1][0] - self.gap
			trace_mat[i][0] = "C"
		for j in range(1,ncols):
			score_mat[0][j] = score_mat[0][j-1] - self.gap
			trace_mat[0][j] = "E"

		for i,x1 in enumerate(seq1):
			for j,x2 in enumerate(seq2):
				values = [score_mat[i][j] + self.__calc_score(x1,x2), score_mat[i+1][j] - self.gap, score_mat[i][j+1] - self.gap]
				direc = "DEC"
				score_mat[i+1][j+1] = max(values)
				trace_mat[i+1][j+1] = direc[values.index(max(values))]

		return score_mat,trace_mat


	def __align(self, seq1: str, seq2: str) -> tuple:
		"""Recebe 2 sequências e retorna o alinhamento das mesmas (Needleman–Wunsch)."""
		nrows = len(seq1) + 1
		ncols = len(seq2) + 1
		score_mat,trace_mat = self.__build_mats(seq1, seq2)
		seq1_aligned = ""
		seq2_aligned = ""
		align_score = score_mat[nrows-1][ncols-1]
		i = nrows - 1
		j = ncols - 1
		while i >= 0 and j >= 0:
			if trace_mat[i][j] == "E":
				seq1_aligned = "-" + seq1_aligned
				seq2_aligned = seq2[j-1] + seq2_aligned
				j -= 1
			elif trace_mat[i][j] == "C":
				seq1_aligned = seq1[i-1] + seq1_aligned
				seq2_aligned = "-" + seq2_aligned
				i -= 1
			elif trace_mat[i][j] == "D":
				seq1_aligned = seq1[i-1] + seq1_aligned
				seq2_aligned = seq2[j-1] + seq2_aligned
				i -= 1
				j -= 1
			else:
				break
		return align_score,seq1_aligned,seq2_aligned


	def __consensus(self, seq1: str, seq2: str) -> str:
		"""Recebe duas sequências alinhadas e devolve a sequência de consenso."""
		return "".join(x1 if x1 != "-" else x2 for x1,x2 in zip(seq1,seq2))


	def mult_align(self) -> list:
		"""Recebe uma lista de sequências e retorna uma lista com as sequências alinhadas."""
		seq1, seq2, *resto = self.seqs
		score, aligned1, aligned2 = self.__align(seq1,seq2)
		alignment = [aligned1,aligned2]
		consenso = self.__consensus(aligned1,aligned2)
		for seq in resto:
			score, aligned1, aligned2 = self.__align(consenso,seq)
			alignment.append(aligned2)
			consenso = self.__consensus(aligned1,aligned2)
		alignment = []
		for seq in self.seqs:
			score, aligned1, aligned2 = self.__align(consenso,seq)
			alignment.append(aligned2)
		return alignment

