"""Quick audit of solubility.py against literature HSP values."""
from topon.singlechain.solubility import estimate_hsp, compute_ra

# Literature HSP values (MPa^0.5) from Hansen Solubility Parameters: A User's Handbook
# and Van Krevelen Properties of Polymers 4th ed.
tests = {
    "toluene":     ("Cc1ccccc1",         18.0,  1.4,  2.0),
    "methanol":    ("CO",                 15.1, 12.3, 22.3),
    "MEK":         ("CCC(C)=O",          16.0,  9.0,  5.1),
    "iso-octane":  ("CC(C)CC(C)(C)C",    14.3,  0.0,  0.0),
    "IPA":         ("CC(O)C",            15.8,  6.1, 16.4),
    "water":       ("O",                 15.5, 16.0, 42.3),
    "PDMS":        ("[Si](C)(C)O",       15.9,  0.0,  4.7),
    "hexane":      ("CCCCCC",            14.9,  0.0,  0.0),
    "acetone":     ("CC(C)=O",           15.5, 10.4,  7.0),
    "ethanol":     ("CCO",               15.8,  8.8, 19.4),
}

print(f"{'Molecule':<14} {'Calc_dD':>8} {'Lit_dD':>8} {'Calc_dP':>8} {'Lit_dP':>8} {'Calc_dH':>8} {'Lit_dH':>8}  {'dD_err%':>8} {'dP_err':>8} {'dH_err':>8}")
print("-" * 110)

for name, (smi, lit_d, lit_p, lit_h) in tests.items():
    hsp = estimate_hsp(smi, verbose=False)
    dd_err = abs(hsp.delta_d - lit_d) / max(lit_d, 0.01) * 100
    dp_err = abs(hsp.delta_p - lit_p)
    dh_err = abs(hsp.delta_h - lit_h)
    flag = " <<<" if dd_err > 30 or dp_err > 5 or dh_err > 10 else ""
    print(f"{name:<14} {hsp.delta_d:8.2f} {lit_d:8.2f} {hsp.delta_p:8.2f} {lit_p:8.2f} {hsp.delta_h:8.2f} {lit_h:8.2f}  {dd_err:7.1f}% {dp_err:8.2f} {dh_err:8.2f}{flag}")

print()
print("=== Verbose breakdown for toluene ===")
estimate_hsp("Cc1ccccc1", verbose=True)
print()
print("=== Verbose breakdown for iso-octane ===")
estimate_hsp("CC(C)CC(C)(C)C", verbose=True)
print()
print("=== Verbose breakdown for methanol ===")
estimate_hsp("CO", verbose=True)
print()
print("=== Verbose breakdown for PDMS ===")
estimate_hsp("[Si](C)(C)O", verbose=True)
