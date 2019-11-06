def hamming_distance(s1, s2):
    i = 0
    count = 0
    while i < len(s1):
        if s1[i] != s2[i]:
            count += 1
        i += 1
    return count

def find_neighbours(pattern, d):
    neighbours = []
    if d == 0:
        neighbours.append(pattern)
        return neighbours
    if len(pattern) == 1:
        return ['A', 'G', 'C', 'T']
    suffix_neighbours = find_neighbours(pattern[1:], d)
    for suffix in suffix_neighbours:
        if hamming_distance(pattern[1:], suffix) < d:
            for nucleotide in ['A', 'T', 'G', 'C']:
                neighbours.append(nucleotide + suffix)
        else:
            neighbours.append(pattern[0] + suffix)
    return neighbours

def complementary_dna_sequence(text):
    complement_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    result = [complement_dict[text[i]] for i in range(len(text) - 1, -1, -1)]
    str = ""
    return str.join(result)

def frequent_mismatching(genome, hamming, k):
    result = []
    patterns = []
    dict = {}
    for i in range(len(genome) - k + 1):
        cur = genome[i:i + k]
        cur_reverse = complementary_dna_sequence(cur)
        result.extend(find_neighbours(cur, hamming))
        result.extend(find_neighbours(cur_reverse, hamming))
    for item in result:
        if item not in dict:
            dict[item] = 1
        else:
            dict[item] += 1
    maximum = 0
    for key in dict.keys():
        if dict[key] > maximum:
            maximum = dict[key]
    for key in dict.keys():
        if dict[key] == maximum:
            patterns.append(key)
    return patterns

if __name__ == "__main__":
    pattern = ""
    genome = "ATCCGTCTACCTCGTCACCTATCCCTCCTCACAATCATCATCCCTATCCGTCCTCACTACACCTCACCTCACCTATCATCATCCACCTATCCCTCTAATCCACCTCTACTACACACGTCTACCTCACACTAATCCACTACACCTATCCCTATCCACTACACTACTAATCCACACGTCGTCCTCACCTCTACGTCTACCTCCTCCTCGTCCTCACTACACA"
    hamming = 2
    k = 6
    print(len(find_neighbours("TGCAT", 2)))
    # print(*frequent_mismatching(genome, hamming, k), sep= ' ')
    # print(hamming_distance("CAGAAAGGAAGGTCCCCATACACCGACGCACCAGTTTA", "CACGCCGTATGCATAAACGAGCCGCACGAACCAGAGAG"))

