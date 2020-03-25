# coding=utf-8
import os

# n = 10
s = 10000
for n in range(10, 11, 10):
    for p in range(2, 11):
        for cp in range(-20, -9):
            for i in range(50):
                command = "./daggen --dot -jump 2 -s " + str(s) + " -n " + str(n) + " -p " + str(p) + " -c " + str(cp)
                os.system(command)
