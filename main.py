from dynamicMSA import align

sequences = ["ACGTACGT", "ACGTTCGT", "ACGTCCGT"]
names     = ["Dizi1", "Dizi2", "Dizi3"]

result = align(sequences, names=names, visualize=True)
print(result)