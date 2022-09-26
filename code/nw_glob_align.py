
class NeedlemanWunsch:

	"""
	Implementa o algoritmo de Needleman-Wunsch para o alinhamento global de sequências biológicas
	"""

	def __init__(self, s1, s2, gap = -1, mol_type = "dna"):

		"""
		Inicializa uma instância da classe Neeedleman-Wunsch
		"""

		self.s1 = s1
		self.s2 = s2
		self.gap = gap
		self.mol_type = mol_type
		if mol_type == "dna": self.scoring = lambda x1, x2: 2 if x1 == x2 else 0
		elif mol_type == "protein":
			blosum_dic = NeedlemanWunsch.__get_blosum_dic("blosum62.txt")
			self.scoring = lambda x1, x2: blosum_dic[x1][x2]
		self.trace_mat, self.score_mat = self.__get_mats()


	@staticmethod
	def __get_blosum_dic(file):

		"""
		Retorna um dicionário de dicionários representativo da matriz blosum62
		"""

		with open(file) as blosum:
			ignore = blosum.readline()
			letras = blosum.readline().split()
			nums = [line.split()[1:] for line in blosum.readlines()]

		blosum_dic = {}
		for i, letra1 in enumerate(letras):
			blosum_dic[letra1] = {}
			for letra2, num in zip(letras, nums[i]):
				blosum_dic[letra1][letra2] = int(num)

		return blosum_dic


	def __get_mats(self):

		"""
		Retorna as matrizes de score e trace construídas a partir das sequências
		"""

		m, n = len(self.s1), len(self.s2)
		trace_mat = [["" for j in range(n + 1)] for i in range(m + 1)]
		score_mat = [[0 for j in range(n + 1)] for i in range(m + 1)]

		for i in range(1, m + 1):
			trace_mat[i][0] = "C"
			score_mat[i][0] = self.gap * i

		for j in range(1, n + 1):
			trace_mat[0][j] = "E"
			score_mat[0][j] = self.gap * i

		direc = "ECD"
		for i in range(1, m + 1):
			for j in range(1, n + 1):
				temp = [
						score_mat[i][j-1] + self.gap, 
						score_mat[i-1][j] + self.gap, 
						score_mat[i-1][j-1] + self.scoring(self.s1[i-1], self.s2[j-1])
					   ]
				score_mat[i][j] = max(temp)
				indices = [i for i in range(len(temp)) if temp[i] == max(temp)]
				trace_mat[i][j] = "".join([direc[i] for i in indices])

		return trace_mat, score_mat


	def __get_unproc_aligns(self, i, j):

		"""
		Retorna uma lista contendo os alinhamentos por processar
		"""

		# base case
		if self.trace_mat[i][j] == "": return [[]]

		# recursive calls
		problem = []
		for p in self.trace_mat[i][j]:
			if p == "E":
				x1, x2 = "-", self.s2[j-1]
				sub_problem = self.__get_unproc_aligns(i, j - 1)
			elif p == "C":
				x1, x2 = self.s1[i-1], "-"
				sub_problem = self.__get_unproc_aligns(i - 1, j)
			elif p == "D":
				x1, x2 = self.s1[i-1], self.s2[j-1]
				sub_problem = self.__get_unproc_aligns(i - 1, j - 1)
			problem += [path + [x1 + x2] for path in sub_problem]

		return problem


	def get_alignments(self):

		"""
		Retorna os melhores alinhamentos de sequências e o respetivo score
		"""

		m, n = len(self.s1), len(self.s2)
		unproc_aligns = self.__get_unproc_aligns(m, n)

		alignments = []
		for align in unproc_aligns:
			s1 = "".join([c[0] for c in align])
			s2 = "".join([c[1] for c in align])
			alignments += [[s1, s2]]

		return alignments, self.score_mat[-1][-1]



if __name__ == "__main__":

	nw_dna = NeedlemanWunsch("ATGAAGGT", "AGAGAGGC", mol_type = "dna")
	print(nw_dna.get_alignments())

	nw_protein = NeedlemanWunsch("GKYESVI", "KYVSSWI", mol_type = "protein")
	print(nw_protein.get_alignments())

