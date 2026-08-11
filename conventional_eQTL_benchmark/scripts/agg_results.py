import pandas as pd

ins = snakemake.input
output = snakemake.output[0]

try:
	results = pd.read_csv(output, sep = "\t")
except FileNotFoundError:
	cols = ['effect_size', 'sample_size', 'cell_count', 'rep', 
         'FN', 'FP', 'TN', 'TP', 'beta_corr']

	results = pd.DataFrame(columns = cols)
	results.to_csv('results/summary_results.tsv',
		sep = "\t", index = None)

for item in ins:
	try:
		row = pd.read_csv(item, sep = "\t")
		if len(results) == 0:
			results = row
		else:
			results = pd.concat([row, results])
	except pd.errors.EmptyDataError:
		print(item)

results.to_csv(output, sep = "\t", index = None)