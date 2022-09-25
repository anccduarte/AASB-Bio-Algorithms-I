
# -*- coding: utf-8 -*-

class Blast:
	"""Implementação de uma versão simplificada do algoritmo de BLAST."""

	def __init__(self, query: str, seq: str, w: int) -> None:

		if type(query) != str or type(seq) != str:
			raise TypeError("A query e a sequência devem ser do tipo 'string'.")

		for c in query.upper() + seq.upper():
			if c not in "ACTG":
				raise ValueError("A query e a sequência devem corresponder a DNA.")

		if type(w) != int:
			raise TypeError("O parâmetro 'w' deve ser do tipo 'int'.")

		if w < 2 or w > len(query):
			raise ValueError("O valor do parâmetro 'w' deve situar-se no intervalo [2, len(query)].")

		self.query = query.upper()
		self.seq = seq.upper()
		self.w = w


	def __str__(self) -> str:
		"""Devolve a query, a sequência e o melhor hit de uma forma mais legível."""
		return f"Query: '{self.query}'\nSequência: '{self.seq}'\nHits: {self.hits()}\nBest hit: {self.best_hit()}"


	def __query_map(self) -> dict:
		"""Recebe a sequência de query e o 'w', e devolve um dicionário cujas keys são sequências e os values são uma lista dos índices."""
		subseq = [self.query[i:i+self.w] for i in range(len(self.query)) if i+self.w <= len(self.query)]
		dic = {}
		for i,s in enumerate(subseq):
			if dic.get(s) == None:
				dic[s] = [i]
			else:
				dic[s].append(i)
		return dic


	def hits(self) -> list:
		"""Recebe o dicionário criado no método '__query_map' e devolve uma lista de hits cujos elementos são tuplos de índices."""
		subseq = [self.seq[i:i+self.w] for i in range(len(self.seq)) if i+self.w <= len(self.seq)]
		qm = self.__query_map()
		hits_list = []
		for key in qm:
			for item in qm[key]:
				for i in range(len(subseq)):
					if key == subseq[i]:
						hits_list.append((item,i))
		return hits_list


	def __extend_hit(self, hit: tuple) -> tuple:
		"""Recebe a query, a sequência, o 'w' e um hit, e devolve um tuplo com o início, o tamanho e o número de matches do hit."""

		query_start = hit[0]
		seq_start = hit[1]
		mismatches = 0
		matches = 1
		extended = 1
		q_inicio_hit = 0
		s_inicio_hit = 0
		go_forward = True
		go_backwards = True
		i = 1

		while go_forward or go_backwards:
			# avançar na query e na sequência
			if query_start+i < len(self.query) and seq_start+i < len(self.seq) and mismatches <= 0.5*extended:
				if self.query[query_start+i] != self.seq[seq_start+i]:
					mismatches += 1
				else:
					matches += 1
				extended += 1
			else:
				go_forward = False
			# retroceder na query e na sequência
			if query_start-i >= 0 and seq_start-i >= 0 and mismatches <= 0.5*extended:
				if self.query[query_start-i] != self.seq[seq_start-i]:
					mismatches += 1
				else:
					matches += 1
				extended += 1
				q_inicio_hit = query_start-i
				s_inicio_hit = seq_start-i
			else:
				go_backwards = False
			i += 1

		return q_inicio_hit,s_inicio_hit,extended,matches


	def __sort_hits(self, lista: list) -> list:
		"""Recebe uma lista de hits e ordena-os por ordem decrescente de score e por ordem crescente de tamanho."""
		return sorted(lista, key = lambda x: (-x[3],x[2]))


	def best_hit(self) -> tuple:
		"""Recebe a query, a sequência e o 'w', e devolve o hit cuja extensão obteve maior score."""
		q_map = self.__query_map()
		hits_list = self.hits()
		if hits_list == []:
			return "Não ocorreu qualquer hit. O valor de 'w' deverá ser ajustado."
		else:
			extended_list = []
			for hit in hits_list:
				extended_list.append(self.__extend_hit(hit))
				extended_sorted = self.__sort_hits(extended_list)
			return extended_sorted[0]

