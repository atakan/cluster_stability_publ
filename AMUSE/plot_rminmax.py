import matplotlib.pyplot as plt

try:
    star_number=int(raw_input('Input(0 to 9):'))
except ValueError:
    print "Not a number"

fo = open('rminmax.txt', 'r')  ## data file

datas = fo.readlines()

seconds = [] ## time
star = [] ## radiis of star

i = 0
t = 0
while(i + 11 < len(datas)):
	
	seconds.append(float(datas[i]))
	#t += 0.004  # will be determined by main code
	star.append(float(datas[i + star_number]))
	i += 11
	

lines = plt.plot(seconds,star)
plt.setp(lines, color='r', linewidth=0.9)
#plt.yscale('log')
plt.title("Radii of Star vs. Time")
plt.xlabel('T')
plt.ylabel('R')
#plt.axis([0, 2000, 0, 3]) ## graph interval
plt.grid(True)
plt.show()

fo.close();
