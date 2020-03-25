# coding=utf-8
import os
#r = 15 # 容量
s = 5000  # 模拟次数
#p = 3  #处理器数量
j = 2  #跳转次数
cp = -15
#n = 10
#cc = 1.5
#for n in [20, 40, 60, 80, 100]:
for n in range(20, 201, 20):
	for p in range(3, 7):
		for cc in [0.1, 0.5, 1.0, 5.0, 10]:
			for i in range(100):
				command = "./daggen  " + " -n " + str(n) + " -p " + str(p) + " -s " + str(s) + " -c " + str(cp) + " --dot " +  " --jump " + str(j) + " --ccr " + str(cc)
				os.system(command)
					#print (command)


