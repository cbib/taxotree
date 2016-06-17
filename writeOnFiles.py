import numpy as np

def writeText(filename,data):
    fo = open("files/" + filename,"w")
    fo.write(data)
    fo.close()

def writeArray(filename,data,header):
    np.savetxt("files/" + filename,data,"%s"," | "," \n",header + " \n","","")

def writeFile(data,header,typeData="text"):
    filename = raw_input("In which file do you want to write it? [Be careful not to choose a name that already exists as it would truncate the existing file]\n")
    if (typeData == "text"):
        writeText(filename,data)
    elif (typeData == "array"):
       	writeArray(filename,data,header)
    else:
        print "Unknown type of data.\n"
    print "-- End of writing\n"

