import periodictable as pt

#Set up initial population

#Read in monomers from input file
readFile = open('smilesText.txt','r')

#create list of monomer SMILES strings
smilesList = []
for line in readFile:
    smilesList = line.split()

readFile.close()

print(smilesList)

#Calculate MW for each monomers
mwList = []

for monomer in smilesList

    #monIter = iter(monomer)
    #if monIter.next().isupper()

    for charcter in monomer
    tempElement = 


#Create inital polymer population


#Loop
    #Selection

    #Mutation

    #Crossover

    #Check

#Print out SMILES strings meeting MW criteria
