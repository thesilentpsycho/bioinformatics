def assignment1(text, pattern):
	count = 0
	for i in range(0, len(text) - len(pattern) + 1):
		if text[i:i + len(pattern)] == pattern:
			count += 1
	return count


def assignment2(text, k):
	dict = {}
	max = 0
	max_list = []
	for i in range(0, len(text) - k + 1):
		temp = text[i:i + k]
		if text[i:i + k] not in dict:
			dict[text[i:i + k]] = 0
		else:
			dict[text[i:i + k]] += 1
		if dict[text[i:i + k]] == max:
			max_list.append(text[i:i + k])
		if (dict[text[i:i + k]] > max):
			max_list.clear()
			max = dict[text[i:i + k]]
			max_list.append(text[i:i + k])
	return max_list

if __name__ == "__main__":
	print(assignment2("CGGAGGACTCTAGGTAACGCTTATCAGGTCCATAGGACATTCA", 3))
