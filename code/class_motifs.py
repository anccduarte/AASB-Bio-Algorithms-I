
# -*- coding: utf-8 -*-

import math

class Motifs:
	"""Determinação de motifs e perfis probabilísticos (PWM e PSSM)."""

	def __init__(self, alignment: list, seq: str, profile: str, pseudocount: int) -> None:

		if type(alignment) != list:
			raise TypeError("O alinhamento deve ser uma lista de strings.")

		if type(seq) != str:
			raise TypeError("A sequência deve ser uma string.")

		if type(profile) != str:
			raise TypeError("O parâmetro 'profile' deve ser do tipo 'str'.")

		if type(pseudocount) != int and type(pseudocount) != float:
			raise TypeError("O parâmetro 'pseudocount' deverá ser do tipo 'int' ou 'float'.")

		for item in alignment:
			if type(item) != str:
				raise TypeError("Os elementos da lista 'alignment' devem ser do tipo 'str'.") 

		for item in alignment:
			for c in item:
				if c.upper() not in "ACTG":
					raise ValueError("Pelo menos uma sequência do alinhamento não corresponde a DNA.")

		for c in seq:
			if c.upper() not in "ACTG":
				raise ValueError("A sequência inserida não corresponde a DNA.")

		if profile not in ["pwm","pssm"]:
			raise ValueError("O parâmetro 'profile' apenas toma os valores 'pwm' ou 'pssm'.")

		self.alignment = alignment
		self.seq = seq.upper()
		self.pseudocount = pseudocount
		self.profile = profile


	def __str__(self) -> str:
		"""Devolve a sequência, o perfil (pwm / pssm) impresso de uma forma legível, e o melhor hit."""
		return (f"Sequência: '{self.seq}'\nPerfil:\n{self.__print_profile(self.__calc_profile())}\n"
				f"Best result: {self.seq_most()[0]} com probabilidade/score de {self.seq_most()[1]:.4f}")


	def __print_profile(self, perfil: list) -> str:
		"""Imprime um perfil (pwm / pssm) de forma legível."""
		string = ""
		bases = sorted(perfil[0].keys())
		tab = [[f"{p[b]:-5.2f}" for b in bases] for p in perfil]
		for p in zip(*([bases] + tab)):
			string += " ".join(p) + "\n"
		return string[:-1]


	def __calc_profile(self) -> list:
		"""Recebe um alinhamento e uma pseudocontagem, e devolve um perfil (pwm / pssm) na forma de lista de dicionários."""
		pwm = []
		total = len(self.alignment) + 4 * self.pseudocount
		for row in zip(*self.alignment):
			row = "".join(row).upper()
			dic = {}
			for base in "ATCG":
				dic[base] = (row.count(base) + self.pseudocount) / total
			pwm.append(dic)
		if self.profile == "pssm":
			for dic in pwm:
				for k in dic:
					dic[k] = math.log(dic[k] * 4, 2)
			pssm = pwm
			return pssm
		return pwm


	def __seq_result(self, subseq: str) -> int:
		"""Recebe uma sequência e um perfil, e retorna a probabilidade (pwm) / score da sequência (pssm)."""
		prfl = self.__calc_profile()
		if self.profile == "pwm":
			prob = 1
			for i in range(len(subseq)):
				prob *= prfl[i][subseq[i]]
			return prob
		else:
			score = 0
			for i in range(len(subseq)):
				score += prfl[i][subseq[i]]
			return score


	def seq_most(self) -> str:
		"""Recebe uma sequência e um perfil, e retorna a subsequência mais provável (pwm) / com maior score (pssm)."""
		prfl = self.__calc_profile()
		seqs = [self.seq[i:i+len(prfl)] for i in range(len(self.seq)) if i+len(prfl) <= len(self.seq)]
		max_result = 0
		max_i = 0
		for i,s in enumerate(seqs):
			result = self.__seq_result(s)
			if result > max_result:
				max_result = result
				max_i = i 
		return seqs[max_i],max_result

