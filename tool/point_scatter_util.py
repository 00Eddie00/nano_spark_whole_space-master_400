def point_scatter(start, end, base, k=1, positive=True):
    result = [start]
    flag = 1 if positive else -1
    i = 0
    while True:
        step = (k * i + base) * flag
        value = result[i] + step
        if positive is True and value > end:
            result[-1] = end
            break
        elif positive is False and value < end:
            result[-1] = end
            break
        else:
            result.append(value)
            i = i + 1
    return result
