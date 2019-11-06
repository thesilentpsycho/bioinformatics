import random
import time
import math
from c1w3 import greedyMotifSearch


def RandomizedMotifSearch(dnas, k):
	theBestMotifs = []
	theBestScore = math.inf
	for i in range(0, 1000):
		motifs = ['GTC', 'CCC', 'ATA', 'GCT']  # [getRandomSubstring(dna, k) for dna in dnas]
		bestMotifs = motifs
		bestScore = greedyMotifSearch.get_score(greedyMotifSearch.build_profile_matrix_with_pseudocounts(bestMotifs))
		while True:
			profile = greedyMotifSearch.build_profile_matrix_with_pseudocounts(motifs)
			motifs = generateMotifsFromProfileAndDnas(profile, dnas, k)
			currScore = greedyMotifSearch.get_score(greedyMotifSearch.build_profile_matrix_with_pseudocounts(motifs))
			if currScore < bestScore:
				bestMotifs = motifs
				bestScore = currScore
			else:
				if bestScore < theBestScore:
					theBestMotifs = bestMotifs
					theBestScore = bestScore
				break
	return theBestMotifs


def generateMotifsFromProfileAndDnas(profile, dnas, k):
	return [greedyMotifSearch.Profile_most_probable_kmer(dna, k, profile) for dna in dnas]


def getRandomSubstring(dna, k):
	random.seed(int(time.time_ns()))
	index = random.randrange(0, len(dna) - k + 1)
	return dna[index:index + k]


if __name__ == "__main__":
	l = """ATGAGGTC
GCCCTAGA
AAATAGAT
TTGTGCTA"""
	k = 3
	d = 2
	t = 20
	dnas = l.strip().split('\n')
	print(*RandomizedMotifSearch(dnas, k), sep='\n')
