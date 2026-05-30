from pathlib import Path

from bioalgs.translation import convolution, convolutionCyclopeptideSequencing, top_m_freq

p = Path(__file__).parent / "datasets" / "real_spectrum.txt" 

spec = list(map(float, p.read_text().split()))

print(
    " ".join(
    map ( str,
	convolutionCyclopeptideSequencing(
		spec,
		40,
		5000,
		False,
		bucket=True,
        bucket_width=4,
		snap_to_aa_masses=True,
		score_tolerance=0.5,
		parent_mass_tolerance=0.5,
	)
    )
    )
)

