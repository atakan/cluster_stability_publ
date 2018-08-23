import matplotlib.pyplot as plt

fo = open('radiis.txt', 'r')  ## data file

datas = fo.readlines()

seconds = [] ## time
one = [] ## shells
two = []
three = []
four = []
five = []
six = []
seven = []
eight = []
nine = []


i = 0
t = 0
while(i < len(datas)):
	
	seconds.append(t)
	t += 0.0004  # will be determined by main code
	
	one.append(float(datas[i]))
	two.append(float(datas[i+1]))
	three.append(float(datas[i+2]))
	four.append(float(datas[i+3]))
	five.append(float(datas[i+4]))
	six.append(float(datas[i+5]))
	seven.append(float(datas[i+6]))
	eight.append(float(datas[i+7]))
	nine.append(float(datas[i+8]))
	
	i += 9
	

lines = plt.plot(seconds,one,seconds,two,seconds,three,seconds,four,seconds,five,seconds,six,seconds,seven,seconds,eight,seconds,nine)
plt.setp(lines, color='b', linewidth=0.6)
#plt.yscale('log')
plt.title("Radii of Shells vs. Time")
plt.xlabel('T')
plt.ylabel('R')
plt.axis([0, 7, 0, 1.5]) ## graph interval
plt.grid(True)
plt.show()

fo.close();
