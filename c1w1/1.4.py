def frequency(text, k, t):
	dict = {}
	result = []
	for i in range(0, len(text) - k + 1):
		temp = text[i:i + k]
		if text[i:i + k] not in dict:
			dict[text[i:i + k]] = 1
		else:
			dict[text[i:i + k]] += 1
	for key in dict:
		if dict[key] >= t:
			result.append(key)
	return result

def clump_problem(genome, k, l, t):
	result = []
	for i in range(0, len(genome) - l + 1):
		old_result = result
		cur_window = genome[i: i + l]
		old_result = frequency(cur_window, k, t)
		result = list(dict.fromkeys(result))
	return result


if __name__ == "__main__":
	k, l, t = 9, 500, 3
	file = open("e_coli.txt", "r")
	genome = file.read().strip()
	file.close()
	print(len(clump_problem(genome, k, l, t)))
