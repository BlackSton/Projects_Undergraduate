number = 1
Max = 0
while number < 1000000:
    chain = 1
    n = number
    while n != 1:
        if n % 2 == 0:
            n = n / 2
            chain += 1
        else:
            n = 3*n + 1
            chain += 1
    if Max < chain:
        Max = chain
        print(number,Max)
    number += 1