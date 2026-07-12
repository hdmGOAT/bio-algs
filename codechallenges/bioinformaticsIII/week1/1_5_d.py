'''
1. If money = 109, DPChange requires a huge array of size 109. Modify the DPChange algorithm so that the array size required does not exceed the value of the largest coin denomination.

2. Recall that our original goal was to make change, not just compute MinNumCoins(money). Modify DPChange so that it not only computes the minimum number of coins but also returns these coins.
'''

from typing import Dict

def DPChange(money: int, coins: list[int]) -> tuple[int, list[int]]:
    largest = max(coins)

    dp: Dict[int, tuple[int, list[int]]] = {
        0: (0, [])
    }

    for amount in range(1, money + 1):
        best_count: int | None = None
        best_solution: list[int] | None = None

        for coin in coins:
            if coin <= amount and (amount - coin) in dp:
                prev_count, prev_solution = dp[amount - coin]
                candidate = prev_count + 1

                if best_count is None or candidate < best_count:
                    best_count = candidate
                    best_solution = prev_solution + [coin]

        if best_count is None or best_solution is None:
            raise ValueError(f"Cannot make change for {amount}.")

        dp[amount] = (best_count, best_solution)

        obsolete = amount - largest
        if obsolete in dp:
            del dp[obsolete]

    return dp[money]
