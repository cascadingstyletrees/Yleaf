import pandas as pd
import numpy as np
import time

def original_logic(pileupfile):
    pileupfile = pileupfile.copy()
    pileupfile.columns = ["chr", "pos", "refbase", "altbase", "reads"]
    pileupfile["pos"] = pileupfile["pos"].astype(int)

    pileupfile["altbase"] = pileupfile["altbase"].str.split(",")
    pileupfile["reads"] = pileupfile["reads"].str.split(",")

    pileupfile["ref_reads"] = pileupfile["reads"].apply(lambda x: x[0])
    pileupfile["alt_reads"] = pileupfile["reads"].apply(lambda x: x[1:])

    pileupfile["alt_reads_dict"] = pileupfile.apply(
        lambda row: dict(zip(row["altbase"], row["alt_reads"])), axis=1
    )
    pileupfile["alt_reads_dict"] = pileupfile["alt_reads_dict"].apply(
        lambda x: {k: int(v) for k, v in x.items()}
    )
    pileupfile["highest_alt_reads"] = pileupfile["alt_reads_dict"].apply(
        lambda x: max(x.values()) if len(x) > 0 else 0
    )
    pileupfile["highest_alt_reads_base"] = pileupfile["alt_reads_dict"].apply(
        lambda x: max(x, key=x.get) if len(x) > 0 else "NA"
    )
    pileupfile["total_reads"] = pileupfile.apply(
        lambda row: int(row["ref_reads"]) + row["highest_alt_reads"], axis=1
    )
    pileupfile["called_ref_perc"] = pileupfile.apply(
        lambda row: round((int(row["ref_reads"]) / row["total_reads"]) * 100, 1)
        if row["total_reads"] > 0
        else 0,
        axis=1,
    )
    pileupfile["called_alt_perc"] = pileupfile.apply(
        lambda row: round((row["highest_alt_reads"] / row["total_reads"]) * 100, 1)
        if row["total_reads"] > 0
        else 0,
        axis=1,
    )

    pileupfile["called_base"] = pileupfile.apply(
        lambda row: row["refbase"]
        if row["called_ref_perc"] >= row["called_alt_perc"]
        else row["highest_alt_reads_base"],
        axis=1,
    )
    return pileupfile

def optimized_logic(pileupfile):
    pileupfile = pileupfile.copy()
    pileupfile.columns = ["chr", "pos", "refbase", "altbase", "reads"]
    pileupfile["pos"] = pileupfile["pos"].astype(int)

    # Convert reads (e.g. "10,5,2") to lists of integers
    # [10, 5, 2]
    reads_list = pileupfile["reads"].str.split(",").apply(lambda x: [int(i) for i in x])

    pileupfile["ref_reads"] = reads_list.str[0]
    alt_reads_list = reads_list.str[1:]

    alt_bases_list = pileupfile["altbase"].str.split(",")

    # To find highest alt read and its base:
    # We can use list comprehensions or vectorized approaches if possible.
    # Since we need both the max value and its index (to get the base),
    # let's see.

    def get_max_alt(bases, reads):
        if not reads:
            return 0, "NA"
        max_val = -1
        max_base = "NA"
        for b, r in zip(bases, reads):
            if r > max_val:
                max_val = r
                max_base = b
        return max_val, max_base

    # Using zip and list comprehension might be faster than row-wise apply
    res = [get_max_alt(b, r) for b, r in zip(alt_bases_list, alt_reads_list)]
    pileupfile["highest_alt_reads"] = [r[0] for r in res]
    pileupfile["highest_alt_reads_base"] = [r[1] for r in res]

    pileupfile["total_reads"] = pileupfile["ref_reads"] + pileupfile["highest_alt_reads"]

    # Vectorized percentage calculation
    # Handle division by zero using where or fillna
    mask = pileupfile["total_reads"] > 0
    pileupfile["called_ref_perc"] = 0.0
    pileupfile.loc[mask, "called_ref_perc"] = (pileupfile.loc[mask, "ref_reads"] / pileupfile.loc[mask, "total_reads"] * 100).round(1)

    pileupfile["called_alt_perc"] = 0.0
    pileupfile.loc[mask, "called_alt_perc"] = (pileupfile.loc[mask, "highest_alt_reads"] / pileupfile.loc[mask, "total_reads"] * 100).round(1)

    pileupfile["called_base"] = pileupfile["highest_alt_reads_base"]
    ref_mask = pileupfile["called_ref_perc"] >= pileupfile["called_alt_perc"]
    pileupfile.loc[ref_mask, "called_base"] = pileupfile.loc[ref_mask, "refbase"]

    return pileupfile

# Generate synthetic data
N = 100000
df = pd.DataFrame({
    0: ["Y"] * N,
    1: np.arange(N),
    2: np.random.choice(["A", "C", "G", "T"], N),
    3: [",".join(np.random.choice(["A", "C", "G", "T"], np.random.randint(1, 4), replace=False)) for _ in range(N)],
    4: [",".join(map(str, np.random.randint(0, 100, np.random.randint(2, 5)))) for _ in range(N)]
})

# Warmup
original_logic(df.head(100))
optimized_logic(df.head(100))

start = time.time()
res_orig = original_logic(df)
end = time.time()
print(f"Original logic took: {end - start:.4f}s")

start = time.time()
res_opt = optimized_logic(df)
end = time.time()
print(f"Optimized logic took: {end - start:.4f}s")

# Verify correctness
# Note: called_base might differ if there are ties and the logic handles them differently,
# but here it seems they should match.
# Also original uses apply(lambda x: max(x, key=x.get)) which returns the FIRST max base.
# My get_max_alt also returns the first max base.

# Drop alt_reads_dict from res_orig as res_opt doesn't have it (or we don't care)
res_orig_cmp = res_orig.drop(columns=["alt_reads_dict", "alt_reads"])
# Wait, original logic keeps altbase and reads as lists, let's compare only the columns we care about
cols_to_cmp = ["highest_alt_reads", "highest_alt_reads_base", "total_reads", "called_ref_perc", "called_alt_perc", "called_base"]

for col in cols_to_cmp:
    pd.testing.assert_series_equal(res_orig[col], res_opt[col], obj=col)

print("Verification successful! Results match.")
