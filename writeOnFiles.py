import numpy as np

def writeText(filename,data):
    try:
        fo = open("files/" + filename + ".taxotree","w")
        fo.write(data)
        fo.close()
    except IOError:
        print "\n/!\ ERROR: Could not open the file. Please try another name."
        raise ValueError

def writeArray(filename,data,header):
    np.savetxt("files/" + filename + ".taxotree",data,"%s"," | "," \n",header + " \n","\n\nEND OF FILE ****","")

def writeFile(data,header,typeData="text"):
    filename = raw_input("In which file do you want to write it? [Be careful not to choose a name that already exists as it would truncate the existing file]\n")
    if (typeData == "text"):
        writeText(filename,data)
    elif (typeData == "array"):
       	writeArray(filename,data,header)
    else:
        print "\n/!\ ERROR: Unknown type of data."
        raise ValueError
    print "-- End of writing\n"

