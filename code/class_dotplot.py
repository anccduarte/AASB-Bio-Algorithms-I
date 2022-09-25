
# -*- coding: utf-8 -*-

class DotPlot:
	"""Construção de dotplots tendo em conta valores de 'window size' e 'stringency'."""

	def __init__(self, seq1: str, seq2: str, w = 0, l = 0) -> None:

 		if type(seq1) != str or type(seq2) != str:
 			raise TypeError("As sequências devem ser do tipo 'string'.")

 		for c in seq1.upper() + seq2.upper():
 			if c.upper() not in "ACTG":
 				raise ValueError("Pelo menos uma das sequências inseridas não corresponde a DNA.")

 		if type(w) != int or type(l) != int:
 			raise TypeError(f"Os parâmetros 'w' e 'l' devem ser do tipo 'int'.")

 		if w < 0 or l < 0:
 			raise ValueError(f"Os valores dos parâmetros 'w' e 'l' devem ser maiores ou iguais a 0.")

 		self.seq1 = seq1.upper()
 		self.seq2 = seq2.upper()
 		self.w = w
 		self.l = l
 		self.mat = []


	def __str__(self) -> str:
		"""Imprime o código utilizado aquando da criação da instância."""
		return f'DotPlot(seq1 = "{self.seq1}", seq2 = "{self.seq2}", w = {self.w}, l = {self.l})'


	def __build_mat(self,nrows: int, ncols: int) -> list:
		"""Recebe um dado número de linhas e colunas e retorna a matriz respetiva."""
		return [[" " for i in range(ncols)] for j in range(nrows)]


	def dotplot(self) -> None:
		"""Recebe 2 sequências e retorna uma matriz de pontos."""

		if self.w != 0 or self.l != 0:
			return "ERRO: Para valores de 'w' e 'l' diferentes de 0, por favor, utilize o método 'dotplot_wl'."

		else:
			nrows = len(self.seq1) + 1
			ncols = len(self.seq2) + 1
			self.mat = self.__build_mat(nrows,ncols)

			for j in range(1,ncols):
				self.mat[0][j] = self.seq2[j-1]

			for i in range(1,nrows):
				self.mat[i][0] = self.seq1[i-1]

			for i in range(1,nrows):
				for j in range(1,ncols):
					if self.seq1[i-1] == self.seq2[j-1]:
						self.mat[i][j] = "*"
			
			lines = [" ".join(line) for line in self.mat]
			
			print("\n".join(lines))


	def dotplot_wl(self) -> None:
		"""Recebe 2 sequências, uma window size (!=0) e um valor de stringency (!=0), e retorna uma matriz de pontos."""

		if self.w == 0 or self.l == 0:
			return "ERRO: Para valores de 'w' e 'l' iguais a 0, por favor, utilize o método 'dotplot'."

		else:
			seq1_sep = [self.seq1[i:i+self.w] for i in range(len(self.seq1)) if i+self.w <= len(self.seq1)]
			seq2_sep = [self.seq2[i:i+self.w] for i in range(len(self.seq2)) if i+self.w <= len(self.seq2)]
			nrows = len(seq1_sep) + self.w
			ncols = len(seq2_sep) + self.w
			self.mat = self.__build_mat(nrows,ncols)

			for i in range(self.w):
				for j in range(self.w,ncols):
					self.mat[i][j] = seq2_sep[j-self.w][i-self.w]

			for i in range(self.w,nrows):
				for j in range(self.w):
					self.mat[i][j] = seq1_sep[i-self.w][j-self.w]

			for i in range(self.w,nrows):
				for j in range(self.w,ncols):
					counter = 0
					for k in range(self.w):
						if seq1_sep[i-self.w][k] == seq2_sep[j-self.w][k]:
							counter += 1
					if counter >= self.l:
						self.mat[i][j] = "*"

			lines = [" ".join(line) for line in self.mat]
			
			print("\n".join(lines))

