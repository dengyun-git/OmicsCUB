#!bin/python
def chop(oldfile,newfile):
	f = open(oldfile,"r")
	sequence = ""
	while(True):
		line = f.readline()
		if(line == ""):break

		if(line[0] == '>'):	
			if(sequence != ""):
				ff = open(newfile,"a")
				ff.write(repr(name) + ";" + "\n") 
				ff.write(repr(sequence) + ";" + "\n")
				ff.close
				
			name = line.split(" ")[0]
			sequence = ""

	 	else:
	 		sequence += line.strip()
        
	ff = open(newfile,"a")
	ff.write(repr(name) + ";" + "\n") 
	ff.write(repr(sequence) + ";" + "\n")
	ff.close 
	return	


f1 = open("FungiNameList.csv")
while(True):
	line = f1.readline()
	if(line == ""):break
	speciesName = line.replace("\n","")
	fileName1 = speciesName + ".cds.all.fa"
	fileName2 = speciesName + ".cds.Fg.fa"
	chop(fileName1,fileName2)
