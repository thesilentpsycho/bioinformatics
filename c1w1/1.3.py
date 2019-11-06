def print_complementary_dna_sequence(text, complement_dict):
	dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
	for i in range(len(text) - 1, -1, -1):
		print(complement_dict[text[i]], end = '')


def index_of_pattern(pattern, genome, start, l):
	result = genome.find(pattern, start)
	if(result == -1):
		return l
	l.append(result)
	return index_of_pattern(pattern, genome, result + 1, l)


if __name__ == "__main__":
	pattern = "ATA"
	genome = "GACGATATACGACGATA"
	result = []
	# dict = {'A' : 'T', 'T' : 'A', 'G' : 'C' , 'C': 'G'}
	# print_complementary_dna_sequence("GATTACA",dict)
	print(*index_of_pattern(pattern, genome, 0, result), sep=' ')
