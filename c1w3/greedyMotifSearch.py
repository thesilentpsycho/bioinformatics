from c1w2 import neighbours
import itertools
import math
import numpy


def hamming_distance(s1, s2):
	i = 0
	count = 0
	while i < len(s1):
		if s1[i] != s2[i]:
			count += 1
		i += 1
	return count


def MotifEnumeration(dnas, k, d):
	patterns = []
	list_of_lists_of_neighbours = []
	for dna in dnas:
		all_neighbours = []
		for i in range(0, len(dna) - k + 1):
			all_neighbours.extend(neighbours.find_neighbours(dna[i:i + k], d))
		list_of_lists_of_neighbours.append(all_neighbours)
	result = set(list_of_lists_of_neighbours[0])
	for s in list_of_lists_of_neighbours[1:]:
		result.intersection_update(s)
	return result


def MedianString(dnas, k, allKmers):
	minsum = math.inf
	pattern = ""
	for kmer in allKmers:
		currSum = sum([min([hamming_distance(dna[i:i + k], kmer) for i in range(len(dna) - k + 1)]) for dna in dnas])
		if minsum > currSum:
			minsum = currSum
			pattern = kmer
	return pattern

def AllMedianStrings(dnas, k, allKmers):
	minsum = math.inf
	pattern = []
	for kmer in allKmers:
		currSum = sum([min([hamming_distance(dna[i:i + k], kmer) for i in range(len(dna) - k + 1)]) for dna in dnas])
		if minsum >= currSum:
			minsum = currSum
			pattern.append(kmer)
	return pattern

def calc_prob(profile, kmer):
	rows = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
	prob = 1.0
	index = 0
	for i in kmer:
		prob = prob * profile[rows[i], index]
		index += 1
	return prob


def Profile_most_probable_kmer(text, k, profile):
	maxprob = -math.inf
	kmerfound = ""
	for i in range(len(text) - k + 1):
		kmer = text[i:i + k]
		prob = calc_prob(profile, kmer)
		if prob > maxprob:
			maxprob = prob
			kmerfound = kmer
	return kmerfound


def generate_all_kmers(k):
	result = []
	bases = ['A', 'T', 'G', 'C']
	result = [''.join(p) for p in itertools.product(bases, repeat=k)]
	return result


def build_profile_matrix(motifs):
	rows = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
	a = numpy.zeros(shape=(4, len(motifs[0])), dtype=float)
	for motif in motifs:
		for idx, letter in enumerate(motif):
			a[rows[letter], idx] = a[rows[letter], idx] + (1.0 / len(motifs))
	return a

def build_profile_matrix_with_pseudocounts(motifs):
	rows = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
	a = numpy.ones(shape=(4, len(motifs[0])), dtype=float)
	for motif in motifs:
		for idx, letter in enumerate(motif):
			a[rows[letter], idx] = a[rows[letter], idx] + 1
	a = a / float(len(motifs) + 4)
	return a

def get_score(profile):
	x = profile.sum(axis=0) - profile.max(axis=0)
	return sum(x)

def GreedyMotifSearch(dnas, k, t):
	best_motifs = [dna[0:k] for dna in dnas]
	best_motifs_matrix = build_profile_matrix(best_motifs)
	y = get_score(best_motifs_matrix)
	first_dna = dnas[0]
	for i in range(len(first_dna) - k + 1):
		kmer = first_dna[i:i+k]
		motifs = []
		motifs.append(kmer)		#1st motif
		for dna in dnas[1:]:
			currProfile = build_profile_matrix(motifs)
			motifs.append(Profile_most_probable_kmer(dna, k, currProfile))
		currScore = get_score(build_profile_matrix(motifs))
		if currScore < get_score(best_motifs_matrix):
			best_motifs = motifs
			best_motifs_matrix = build_profile_matrix(best_motifs)
	return best_motifs

def GreedyMotifSearchWithPseudocounts(dnas, k, t):
	best_motifs = [dna[0:k] for dna in dnas]
	best_motifs_matrix = build_profile_matrix_with_pseudocounts(best_motifs)
	y = get_score(best_motifs_matrix)
	first_dna = dnas[0]
	for i in range(len(first_dna) - k + 1):
		kmer = first_dna[i:i+k]
		motifs = []
		motifs.append(kmer)		#1st motif
		for dna in dnas[1:]:
			currProfile = build_profile_matrix_with_pseudocounts(motifs)
			motifs.append(Profile_most_probable_kmer(dna, k, currProfile))
		currScore = get_score(build_profile_matrix_with_pseudocounts(motifs))
		if currScore < get_score(best_motifs_matrix):
			best_motifs = motifs
			best_motifs_matrix = build_profile_matrix_with_pseudocounts(best_motifs)
	return best_motifs

if __name__ == "__main__":
	l = """CTCGATGAGTAGGAAAGTAGTTTCACTGGGCGAACCACCCCGGCGCTAATCCTAGTGCCC
GCAATCCTACCCGAGGCCACATATCAGTAGGAACTAGAACCACCACGGGTGGCTAGTTTC
GGTGTTGAACCACGGGGTTAGTTTCATCTATTGTAGGAATCGGCTTCAAATCCTACACAG"""
	k = 7
	d = 2
	t = 25
	dnas = l.strip().split('\n')
	# print(*MotifEnumeration(dnas, k , d), sep= ' ')
	print(*AllMedianStrings(dnas, k, generate_all_kmers(k)))
	temp = """0.276 0.263 0.276 0.224 0.211 0.276 0.276 0.263 0.224 0.276 0.316 0.263 0.237
0.263 0.224 0.224 0.25 0.25 0.237 0.316 0.303 0.276 0.303 0.237 0.289 0.316
0.263 0.224 0.237 0.158 0.171 0.184 0.184 0.263 0.237 0.25 0.211 0.211 0.237
0.197 0.289 0.263 0.368 0.368 0.303 0.224 0.171 0.263 0.171 0.237 0.237 0.211
"""
	profile = numpy.array([x.strip().split(' ') for x in temp.strip().split('\n')]).astype(numpy.float)
	# text = "AAGCGCTCCTATCAACGTAGGTGTATACTCGCTATCAGTTATATGAAGACAATGCCCGGTATTCCATATATATTAGCGACCTCCGTGTTCACCAAGACTACCTTCCGCTTACAAACTATGGACTGTAATTTCTTCGTGATGCGTGTAGCCGTTCCGGGCGGGTCGAACCCTAAGGAGAGAACTCATGAATACCACATAGTGCAGGGTTAAGTGGATTGAATACTTCAACTTCACCGCATTGTGTACCCGACGGGGGGGGTATTGACGAGGGACCATGAATACCTCGTACTGTGCGGAGATAGGAATCTATGAGCCACATCGGCCTGGGAGCTGAACAAAGGCGGTCAAACGTCCGATTATCGGCCATCCGTCTCATGCCCCGGTTCAAATGTTTGTCTGGTAAGCCAGACGCGAGGTGGAACTATGTACCTATGCCGTGGCTTCCCGGCCTTGGCACTCGGCAGATTGTCGAGGTCACACTCACAACGCGAATCTTGTGACCTGATTCGAACTATTTGCTTAGCACGTAAGCAATGAATGAATAACAAGCGTGCTTACCCTAACAAACGTTAGAAGTCTAAGGACCACAACAAAATATCAATCACCATGCAGAATTACTGCCTCACAGAAACTTGTAGAGGAGCGGGGAACCTCTATGTGCTGATTCGTTTGGTACAAATTTTGCCCGGGGGTTGATTCGTATTCACTACGAGGTCTGTTTGAGCATAGCCCCTTCGTGCGTTTGACAACGGGTGGTAGGCAGTACACGGCGATGAGCTGTGCCGCTGAAGACTCCGGCATGCTTAAGTCATACTTGACGACGCTTTCTGCCTAACGAAATTGTCCGAATAGACGACGATGAAATCCACATGACGGGGATAACTTCTGTATACGCCCCATTCAGAAGGCTACCATCAAACTTGAGCCTTAGTCCCCACTGTTACTCTGACTAGAGAATGAGGCGGCAGGGTATTCATAGTGTGGTTGATAGACACTCGCG"

	# print(Profile_most_probable_kmer(text, k, profile))
	# print(calc_prob(profile, "CCGAG"))
	# print(*GreedyMotifSearchWithPseudocounts(dnas, k, t), sep='\n')
