from pathlib import Path
def DPChange(m, Coins):
    tab = {0: 0}
    for i in range(m + 1):
        min = float('inf')
        for c in Coins:
            if i - c >= 0:
                now = tab[i-c] + 1
                if now < min:
                    min = now
                    tab[i] = now
    return tab[m]

path = Path(__file__).parent / "datasets" / "dataset_30195_10.txt"

m_str, coins_str = path.read_text().splitlines()

m = int(m_str)
coins = list(map(int, coins_str.split()))

print(DPChange(m, coins))
