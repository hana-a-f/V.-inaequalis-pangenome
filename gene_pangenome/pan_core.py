import pandas as pd
import itertools
import multiprocessing as mp

def process_combinations(chunk, orth_bin):
    """
    Process a chunk of combinations to calculate core genes, pan genes, and genome length.
    """
    results = []
    for combo in chunk:
        sub_matrix = orth_bin[list(combo)]
        core_genes = (sub_matrix.sum(axis=1) == len(combo)).sum()
        pan_genes = (sub_matrix.sum(axis=1) > 0).sum()
        genome_length = len(combo)
        results.append({
            "Combination": ",".join(combo),
            "Core_Genes": core_genes,
            "Pan_Genes": pan_genes,
            "Genome_Length": genome_length
        })
    return results

if __name__ == "__main__":
    # Step 1: Load the data
    orth_bin = pd.read_csv("Orthogroups.GeneCount.tsv", sep="\t", header=0)

    # Step 2: Keep only the desired columns
    orth_bin = orth_bin.iloc[:, 1:19]

    # Step 3: Convert to binary (1 if orthogroup is present, 0 otherwise)
    orth_bin = (orth_bin > 0).astype(int)

    # Step 4: Generate all combinations of the genomes (columns)
    column_ids = orth_bin.columns
    combinations = []
    for i in range(1, len(column_ids) + 1):
        combinations.extend(itertools.combinations(column_ids, i))

    # Step 5: Split combinations into chunks for multiprocessing
    num_cores = mp.cpu_count()  # Automatically detects the number of available cores
    chunk_size = len(combinations) // num_cores + 1
    chunks = [combinations[i:i + chunk_size] for i in range(0, len(combinations), chunk_size)]

    # Step 6: Initialize multiprocessing pool and process combinations
    with mp.Pool(processes=num_cores) as pool:
        results = pool.starmap(process_combinations, [(chunk, orth_bin) for chunk in chunks])

    # Step 7: Flatten results
    results = [item for sublist in results for item in sublist]

    # Step 8: Convert results to a DataFrame
    results_df = pd.DataFrame(results)

    # Step 9: Save the results to a CSV file
    results_df.to_csv("results.csv", index=False)

