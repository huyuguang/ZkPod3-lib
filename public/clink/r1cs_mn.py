import math
import sys

def f(S,c,m):
    return (S + 1.0 * c * S / m) / (math.log(S / m,10))

def min(S,c):
    min_x = 0
    min_y = -1
    loop = int(S / 10) - 1
    for index in range(loop):
        #print (index)
        x = (index + 1) * 10
        y = f(S,c,x)
        if min_y < 0 or y < min_y:
            min_x = x
            min_y = y
    print(S, c, "min_x: ", min_x, "min_y: ", min_y)
    return
    
def main(argv):
    c = int(argv[1])
    S = [10000, 100000, 1000000,10000000, 100000000,1000000000]
    for i in S:
        min(i,c)
    return

if __name__ == "__main__":
    main(sys.argv)
