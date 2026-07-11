def RecursiveChange(money, Coins, memo=None):
    if memo is None:
        memo = {}

    if money == 0:
        return 0

    if money in memo:
        return memo[money]

    minNumCoins = float('inf')
    for i in range(len(Coins)):
        if money >=  Coins[i]:
            numCoins = RecursiveChange(money - Coins[i], Coins, memo)
            if numCoins + 1 < minNumCoins:
                minNumCoins = numCoins + 1
    memo[money] = minNumCoins
    return minNumCoins

print(RecursiveChange(76, [5,4,1]))
