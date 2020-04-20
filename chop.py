import sys

def chop(oldfile,newfile)

f=fopen(oldfile) 
sequence = ""
while(True):
	line = f.readline()
	    if(line == ""): break 
	    if(line[0] == '>'): 
	    	if(sequence != ""): 
			print name + ";"
			print sequence + ";"
			name = line.split(" ")[0] 
			sequence = ""
		else: 
			sequence += line.strip()

ff=open(newfile,'w')

f3.write('name + ";"')

f3.write('sequence + ";"')

return


