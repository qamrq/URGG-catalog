"""Layout-test stand-in for 63_genome_thiamine_guilds.tsv. NOT real data; only used to test the figure layout."""
import numpy as np
import pandas as pd

rng = np.random.default_rng(1)
sh = pd.read_csv(r"D:\QRQ_ZJU\MAG-SAG\Manuscript\Figures\Fig2\supp\representative_host_range.csv",
                 usecols=["representative", "host_species_list"])
sh = sh.rename(columns={"representative": "genome"}).drop_duplicates("genome").head(25451)
n = len(sh)
tiers = np.where(np.arange(n) < 454, "Core / broad (>3 host species)",
                 np.where(np.arange(n) < 454 + 4811, "Shared (2-3 host species)", "Host-specific (1 host species)"))
guilds = ["1_de_novo", "2_partial_de_novo", "3_precursor_salvage", "4_thiamine_salvage", "5_none_detected"]
share = {"Core / broad (>3 host species)": [.12, .19, .14, .33, .22],
         "Shared (2-3 host species)": [.044, .16, .12, .39, .28],
         "Host-specific (1 host species)": [.025, .10, .10, .44, .33]}
guild = np.array([rng.choice(guilds, p=np.array(share[t]) / sum(share[t])) for t in tiers])
thic = np.isin(guild, guilds[:2]).astype(int)
out = pd.DataFrame({
    "genome": sh.genome.values, "guild": guild, "sharing_category": tiers,
    "hmp_denovo": thic, "transport": rng.random(n) < 0.5, "tpk_thiN": rng.random(n) < 0.6,
    "tmp_kinase_thiL": rng.random(n) < 0.3})
out["transport"] = out.transport.astype(int)
out["tpk_thiN"] = out.tpk_thiN.astype(int)
out["tmp_kinase_thiL"] = out.tmp_kinase_thiL.astype(int)
out.to_csv(r"D:\Claude_workzone\test\mock_63_genome_thiamine_guilds.tsv", sep="\t", index=False)
print(len(out), out.guild.value_counts().to_dict())
