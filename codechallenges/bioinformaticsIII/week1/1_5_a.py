def RecursiveChange(money, Coins):
    if money == 0:
        return 0
    minNumCoins = float('inf')
    for i in range(len(Coins)):
        if money >=  Coins[i]:
            numCoins = RecursiveChange(money - Coins[i], Coins)
            if numCoins + 1 < minNumCoins:
                minNumCoins = numCoins + 1
    return minNumCoins

print(RecursiveChange(76, [5,4,1]))
