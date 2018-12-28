#import periodictable as pt

#Set up initial population

#Read in monomers from input file
readFile = open('../input_files/smilesText.txt','r')

#create list of monomer SMILES strings
smilesList = []
for line in readFile:
    smilesList = line.split()

readFile.close()

print(smilesList)
'''
#Calculate MW for each monomers
mwList = []

for monomer in smilesList
    mwList.append()

    #monIter = iter(monomer)
    #if monIter.next().isupper()

    #for charcter in monomer
        #if character.isupper()
            #if tempElement

            #tempElement = character
        #else
            #tempElement = tempElement + character

    for character in monomer
        if character.isUpper()
            mwList[monomer.index()] += pt.character.mass()


#Create inital polymer population


#Loop
    #Selection

    #Mutation

    #Crossover

    #Check

#Print out SMILES strings meeting MW criteria
'''
