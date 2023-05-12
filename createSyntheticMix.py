import time, os, re

# region: Functions

def makeFlatFileDB(dir: str, regex: str = None, fileExt: str = None, outFile: str = None, maxFilesRead:int = 100000000, excludeDirs: list[str] = [], verbose: bool = True):
    """Finds all files that fit a regex in a specified folder
    :param dir: Directory to search
    :param regex: The regex to search by
    :param outFile: The output file path
    :param maxFilesRead: The maximum number of files to read, defaults to 1000
    :param excludeDirs: List of folders to exclude
    :param verbose: Print progress messages?, defaults to True
    """    
    nFasta = nFiles = speed = 0
    out = []
    lastCheck = startTime = time.time()
    excludeDirs = "|".join(excludeDirs)

    if outFile is not None:
        out = open(outFile,'w')

    if (verbose): print("Searching for files...")

    for root, dirs, files in os.walk(dir, topdown=True):
        dirs.sort(reverse=True)
        if nFasta >= maxFilesRead: break
        for file in files:
            nFiles += 1
            match = file.endswith((fileExt)) if (regex is None) else re.match(regex, file)
            if match:
                if re.search(excludeDirs, root) != None: continue
                path = str(root) + "/" + str(file)
                if outFile is not None: 
                    out.write(path + "\n")
                else: 
                    out.append(path)
                nFasta += 1
            if (nFiles % 1000 == 0): 
                speed = str(round(1000/(time.time() - lastCheck)))
                lastCheck = time.time()
            if (verbose): print("   Parsed {} files and found {} FASTA files ({}/s)                ".format(nFiles,nFasta,speed), end="\r")

    if (verbose): print("   Parsed {} files and found {} FASTA files ({}s)                 ".format(nFiles,nFasta,str(round(time.time() - startTime,2))), end="\r")

    return (out if outFile is None else outFile)

# endregion

os.chdir(os.path.dirname(__file__))
path = "/nfs/APL_Genomics/virus_covid19/routine_seq/2023_01_Runs"
BAMfiles = "./BAMfiles.txt"

if (not os.path.isfile(BAMfiles)):
    makeFlatFileDB(path, fileExt=((".bam")), outFile="BAMfiles", excludeDirs=["work","tmp_bam","qc","test","troubleshooting"])

BAMfiles = [line.strip() for line in open(BAMfiles, 'r')]




