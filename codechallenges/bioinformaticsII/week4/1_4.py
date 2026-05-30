from pathlib import Path

from bioalgs.translation import convolution, convolutionCyclopeptideSequencing, top_m_freq

p = Path(__file__).parent / "datasets" / "real_spectrum.txt" 

spec = list(map(float, p.read_text().split()))
spec = [round(x) for x in spec]

print(convolutionCyclopeptideSequencing(spec, 40, 5000, True, bucket=True))
