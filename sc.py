import os

#os.system("./init  1000 1 0")
#os.system("./main")

n = 0.5
while n <= 1.1:
	m = -0.3
	while m <= 0 :
		script_output = open('script_output.txt' , 'a')
		script_output.write("n = " + str(n) + "   m = " + str(m) + "\n")
		script_output.close()
		for num in range(0,10):
			os.system("./init  1000 " + str(n) + "  " + str(m))
			os.system("./main")
		m += 0.1
	n += 0.1