
# -*- coding: utf-8 -*-

class Align:
	"""Construção de alinhamentos locais (Smith-Waterman) e globais (Needleman-Wunsch)."""

	def __init__(self, seq1: str, seq2: str, align_type = "global", scoring = [2,0], gap = 4) -> None:

		if type(seq1) != str or type(seq2) != str:
			raise TypeError("As sequências a alinhar devem ser do tipo 'string'.")

		if align_type not in ["global","local"]:
			raise ValueError(f"Um alinhamento do tipo '{align_type}' não é válido.")

		if ((type(scoring) == list and (len(scoring) != 2 or type(scoring[0]) != int or type(scoring[1]) != int)) or
		   (type(scoring) == str and scoring not in ["blosum50","blosum62","blosum80"]) or 
		   (type(scoring) != list and type(scoring) != str)):
			raise TypeError("O parâmetro 'scoring' apenas aceita uma lista de dois valores inteiros ou uma matriz 'blosum(50|62|80)'.")

		if type(gap) != int:
			raise TypeError("O parâmetro 'gap' deve ser do tipo 'int'.")

		if gap < 0:
			raise ValueError("O valor do parâmetro 'gap' não deve ser menor que 0.")


		self.seq1 = seq1.upper()
		self.seq2 = seq2.upper()
		self.align_type = align_type
		self.scoring = scoring
		self.gap = gap


	def __str__(self) -> str:
		"""Devolve as sequências alinhadas e o score do alinhamento impressos de uma forma mais legível."""
		seq1,seq2,score = self.align()
		return f"Sequence 1 aligned: {seq1}\nSequence 2 aligned: {seq2}\nAlignment score: {score}"


	def __build_dic(self) -> dict:
		"""Converte uma matriz blosum contida num ficheiro txt num dicionário de dicionários."""
		matrix = open(f"{self.scoring}.txt")
		dic = {}
		headers,*mat = [linha.split() for linha in matrix.readlines() if linha.strip() != ""]
		matrix.close()
		for lin in mat:
			letra,*num = lin
			dic[letra] = {}
			for outra_letra,valor in zip(headers,num):
				dic[letra][outra_letra] = valor
		return dic


	def __calc_score(self, x1: str, x2: str) -> int:
		"""Recebe dois caracteres e devolve o score correspondente."""
		if type(self.scoring) == list:
			return self.scoring[0] if x1 == x2 else self.scoring[1]
		else:
			dic_blosum = self.__build_dic()
			return int(dic_blosum[x1][x2])


	def __build_mats(self) -> list:
		"""Recebe 2 sequências, um tipo de alinhamento e o valor da gap penalty, e retorna as matrizes de score e trace."""

		nrows = len(self.seq1) + 1
		ncols = len(self.seq2) + 1
		score_mat = [[0 for j in range(ncols)] for i in range(nrows)]
		trace_mat = [[0 for j in range(ncols)] for i in range(nrows)]

		if self.align_type == "global":
			for i in range(1,nrows):
				score_mat[i][0] = score_mat[i-1][0] - self.gap
			for j in range(1,ncols):
				score_mat[0][j] = score_mat[0][j-1] - self.gap

		for i in range(1,nrows):
			trace_mat[i][0] = "C"
		for j in range(1,ncols):
			trace_mat[0][j] = "E"

		for i,x1 in enumerate(self.seq1):
			for j,x2 in enumerate(self.seq2):
				if self.align_type == "global":
					values = [score_mat[i+1][j] - self.gap, score_mat[i][j+1] - self.gap, score_mat[i][j] + self.__calc_score(x1,x2)]
					direc = "ECD"
				elif self.align_type == "local":
					values = [score_mat[i+1][j] - self.gap, score_mat[i][j+1] - self.gap, score_mat[i][j] + self.__calc_score(x1,x2), 0]
					direc = "ECD0"
				score_mat[i+1][j+1] = max(values)
				trace_mat[i+1][j+1] = direc[values.index(max(values))]

		return score_mat,trace_mat


	def __max_score(self) -> tuple:
		"""Recebe uma matriz de score, e retorna o valor e os índices do score máximo na mesma."""
		score_mat = self.__build_mats()[0]
		score = 0
		i_max = 0
		j_max = 0
		for i in range(len(score_mat)):
			for j in range(len(score_mat[0])):
				if score_mat[i][j] > score:
					score = score_mat[i][j]
					i_max = i
					j_max = j
		return score,i_max,j_max


	def align(self) -> tuple:
		"""Recebe 2 sequências e um tipo de alinhamento, e retorna as sequências alinhadas e o score do alinhamento."""

		nrows = len(self.seq1) + 1
		ncols = len(self.seq2) + 1
		score_mat,trace_mat = self.__build_mats()
		seq1_aligned = ""
		seq2_aligned = ""

		if self.align_type == "global":
			align_score = score_mat[nrows-1][ncols-1]
			i = nrows - 1
			j = ncols - 1
			while i >= 0 and j >= 0:
				if trace_mat[i][j] == "E":
					seq1_aligned = "-" + seq1_aligned
					seq2_aligned = self.seq2[j-1] + seq2_aligned
					j -= 1
				elif trace_mat[i][j] == "C":
					seq1_aligned = self.seq1[i-1] + seq1_aligned
					seq2_aligned = "-" + seq2_aligned
					i -= 1
				elif trace_mat[i][j] == "D":
					seq1_aligned = self.seq1[i-1] + seq1_aligned
					seq2_aligned = self.seq2[j-1] + seq2_aligned
					i -= 1
					j -= 1
				else:
					break

		elif self.align_type == "local":
			align_score,i_start,j_start = self.__max_score()
			while i_start >= 0 and j_start >= 0:
				if score_mat[i_start][j_start] != 0:
					if trace_mat[i_start][j_start] == "E":
						seq1_aligned = "-" + seq1_aligned
						seq2_aligned = self.seq2[j_start-1] + seq2_aligned
						j_start -= 1
					elif trace_mat[i_start][j_start] == "C":
						seq1_aligned = self.seq1[i_start-1] + seq1_aligned
						seq2_aligned = "-" + seq2_aligned
						i_start -= 1
					elif trace_mat[i_start][j_start] == "D":
						seq1_aligned = self.seq1[i_start-1] + seq1_aligned
						seq2_aligned = self.seq2[j_start-1] + seq2_aligned
						i_start -= 1
						j_start -= 1
				else:
					break

		return seq1_aligned,seq2_aligned,align_score

