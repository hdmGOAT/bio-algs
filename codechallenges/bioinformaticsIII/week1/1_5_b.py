def MinNumCoins(m, Coins):
    tab = {0: 0}
    for i in range(m + 1):
        min = float('inf')
        for c in Coins:
            if i - c >= 0:
                now = tab[i-c] + 1
                if now < min:
                    min = now
                    tab[i] = now
    return list(tab.values())
l = MinNumCoins(22, [1,4,5])
print(" ".join(map(str, l[13:])))
