from dynamicMSA.align import align
sequences = ["ACGTACGT", "ACGGTTCGT", "CGTCCGT","ATGCGGTA"]
names     = ["Dizi1", "Dizi2", "Dizi3","Dizi4"]
result = align(sequences, names=names, visualize=True)