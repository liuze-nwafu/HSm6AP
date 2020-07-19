import os
import sys
import re
import itertools
from collections import Counter
import re
import itertools
from collections import Counter
def kmerArray(sequence, k):
    kmer = []
    for i in range(len(sequence) - k + 1):
        kmer.append(sequence[i:i + k])
    return kmer

def train2Kmer(fastas, k=4, type="DNA", upto=False, normalize=True, **kw):
    encoding = []
    header = ['class']
    NA = 'ACGT'
    if type in ("DNA", 'RNA'):
        NA = 'ACGT'
    else:
        NA = 'ACDEFGHIKLMNPQRSTVWY'

    if k < 1:
        print('Error: the k-mer value should larger than 0.')
        return 0

    if upto == True:
        for tmpK in range(1, k + 1):
            for kmer in itertools.product(NA, repeat=tmpK):
                header.append(''.join(kmer))
        encoding.append(header)
        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])
            count = Counter()
            for tmpK in range(1, k + 1):
                kmers = kmerArray(sequence, tmpK)
                count.update(kmers)
                if normalize == True:
                    for key in count:
                        if len(key) == tmpK:
                            count[key] = count[key] / len(kmers)
            code = [1]
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    else:
        for kmer in itertools.product(NA, repeat=k):
            header.append(''.join(kmer))
        encoding.append(header)
        for i in fastas:
            name, sequence = i[0], re.sub('-', '', i[1])
            kmers = kmerArray(sequence, k)
            count = Counter()
            count.update(kmers)
            if normalize == True:
                for key in count:
                    count[key] = count[key] / len(kmers)
            code = [1]
            for j in range(2, len(header)):
                if header[j] in count:
                    code.append(count[header[j]])
                else:
                    code.append(0)
            encoding.append(code)
    return encoding

def train2CKSNAP(fastas, gap=5, **kw):
    if gap < 0:
        print('Error: the gap should be equal or greater than zero' + '\n\n')
        return 0

    AA = 'ACGT'
    encodings = []
    aaPairs = []
    for aa1 in AA:
        for aa2 in AA:
            aaPairs.append(aa1 + aa2)

    header = ['class']
    for g in range(gap + 1):
        for aa in aaPairs:
            header.append(aa + '.gap' + str(g))
    encodings.append(header)

    for i in fastas:
        name, sequence = i[0], i[1]
        code = [1]
        for g in range(gap + 1):
            myDict = {}
            for pair in aaPairs:
                myDict[pair] = 0
            sum = 0
            for index1 in range(len(sequence)):
                index2 = index1 + g + 1
                if index1 < len(sequence) and index2 < len(sequence) and sequence[index1] in AA and sequence[
                    index2] in AA:
                    myDict[sequence[index1] + sequence[index2]] = myDict[sequence[index1] + sequence[index2]] + 1
                    sum = sum + 1
            for pair in aaPairs:
                code.append(myDict[pair] / sum)
        encodings.append(code)
    return encodings

def read_nucleotide_sequences(file):
    if os.path.exists(file) == False:
        print('Error: file %s does not exist.' % file)
        sys.exit(1)
    with open(file) as f:
        records = f.read()
    if re.search('>', records) == None:
        print('Error: the input file %s seems not in FASTA format!' % file)
        sys.exit(1)
    records = records.split('>')[1:]
    fasta_sequences = []
    for fasta in records:
        array = fasta.split('\n')
        header, sequence = array[0].split()[0], re.sub('[^ACGTU-]', '-', ''.join(array[1:]).upper())
        header_array = header.split('|')
        name = header_array[0]
        label = header_array[1] if len(header_array) >= 2 else '0'
        label_train = header_array[2] if len(header_array) >= 3 else 'training'
        sequence = re.sub('U', 'T', sequence)
        fasta_sequences.append([name, sequence, label, label_train])
    return fasta_sequences

import re
def train2DNC(fastas, **kw):
    base = 'ACGT'

    encodings = []
    dinucleotides = [n1 + n2 for n1 in base for n2 in base]
    header = ['class'] + dinucleotides
    encodings.append(header)

    AADict = {}
    for i in range(len(base)):
        AADict[base[i]] = i

    for i in fastas:
        name, sequence = i[0], re.sub('-', '', i[1])
        code = [1]
        tmpCode = [0] * 16
        for j in range(len(sequence) - 2 + 1):
            tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] = tmpCode[AADict[sequence[j]] * 4 + AADict[sequence[j+1]]] +1
        if sum(tmpCode) != 0:
            tmpCode = [i/sum(tmpCode) for i in tmpCode]
        code = code + tmpCode
        encodings.append(code)
    return encodings


chemical_property = {
    'A': [1, 1, 1],
    'C': [0, 1, 0],
    'G': [1, 0, 0],
    'T': [0, 0, 1],
    'U': [0, 0, 1],
    '-': [0, 0, 0],
}


def train2NCP(fastas, **kw):

    AA = 'ACGT'
    encodings = []
    header = ['class']
    for i in range(1, len(fastas[0][1]) * 3 + 1):
        header.append('NCP.F'+str(i))

    for i in fastas:
        name, sequence= i[0], i[1]
        code = [1]
        for aa in sequence:
            code = code + chemical_property.get(aa, [0, 0, 0])
        encodings.append(code)
    return encodings

def main1(name,TestPosName):
    pathTest = name + "Data/"
    pathTestOut = name + "temp/"
    fileTest1 = train2DNC(read_nucleotide_sequences(pathTest + TestPosName))
    fea = "DNC.csv"
    dfTest1 = pd.DataFrame(fileTest1)
    dfTest1.fillna(0)
    dfTest1.to_csv(pathTestOut + TestPosName[:-6] + fea, header=None, index=None)
    print(fea + "end.....")

def main2(name,TestPosName):
    pathTest = name + "Data/"
    pathTestOut = name + "temp/"
    fileTest1 = train2Kmer(read_nucleotide_sequences(pathTest + TestPosName))
    fea = "Kmer.csv"
    dfTest1 = pd.DataFrame(fileTest1)
    dfTest1.fillna(0)
    dfTest1.to_csv(pathTestOut +TestPosName[:-6] + fea, header=None, index=None)
    print(fea + "end.....")

def main3(name,TestPosName):
    pathTest = name + "Data/"
    pathTestOut = name + "temp/"
    fileTest1 = train2CKSNAP(read_nucleotide_sequences(pathTest + TestPosName))
    fea = "CKSNAP.csv"
    dfTest1 = pd.DataFrame(fileTest1)
    dfTest1.fillna(0)
    dfTest1.to_csv(pathTestOut + TestPosName[:-6] + fea, header=None, index=None)
    print(fea + "end.....")

def main4(name,TestPosName):
    pathTest = name + "Data/"
    pathTestOut = name + "temp/"
    fileTest1 = train2NCP(read_nucleotide_sequences(pathTest + TestPosName))
    fea = "NCP.csv"
    dfTest1 = pd.DataFrame(fileTest1)
    dfTest1.fillna(0)
    dfTest1.to_csv(pathTestOut + TestPosName[:-6]+ fea, header=True, index=None)
    print(fea + "end.....")

def orderPse(name,TestPosName):
    os.system(
        "python "+name+"pseInOne/pse.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA PC-PseDNC-General -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"PC-PseDNC.txt -w 0.1 -k 1")
    os.system(
        "python "+name+"pseInOne/pse.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA PC-PseTNC-General -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"PC-PseTNC.txt -w 0.1 -k 1")
    os.system(
        "python "+name+"pseInOne/pse.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA SC-PseTNC-General -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"SC-PseTNC.txt -w 0.1 -k 1")
    os.system(
        "python "+name+"pseInOne/pse.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA SC-PseDNC-General -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"SC-PseDNC.txt -w 0.1 -k 1")
    os.system(
        "python "+name+"pseInOne/pse.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA PseDNC -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"PseDNC.txt -w 0.1 -k 1")
    os.system(
        "python "+name+"pseInOne/nac.py " + name + "model/2.fasta " + name + "Data/"+TestPosName+" DNA Mismatch -k 3 -f tab -out " + name + "temp/2.txt " + name + "temp/"+TestPosName[:-6]+"Mismatch.txt")

import numpy as np
import pandas as pd
def txt2csv(name,TestPosName):
    features=["Mismatch","SC-PseTNC","SC-PseDNC","PC-PseTNC","PC-PseDNC"]
    for strFea in features:
        pathTest=name+"temp/"
        TesttxtPos = np.loadtxt(pathTest + TestPosName[:-6] + strFea + ".txt")
        header = ["class"]
        for i in range(len(TesttxtPos[0])):
            header.append(str(i))
        classCol = np.array(list(np.ones(len(TesttxtPos),int)))
        # TesttxtPos = np.insert(TesttxtPos, 0, values=classCol, axis=1)
        TesttxtPos = np.insert(TesttxtPos, 0, values=classCol, axis=1)
        TesttxtDF1 = pd.DataFrame(TesttxtPos)
        TesttxtDF1.to_csv(name + "temp/"+TestPosName[:-6] + strFea + ".csv", index=False, encoding='utf-8', header=header)
        print(strFea + "test_end.......")

import csv
import numpy as np

def colNum2(file):
    with open(file, "r") as csvfile:
        reader = csv.reader(csvfile)
        rows=[row for row in reader]
        a=[]
        for i in range(len(rows[0])):
            if rows[0][i] =="class":
                continue
            else:
                a.append(rows[0][i])
        return a
def colReal2(file,index):
    with open(file, "r") as csvfile:
        reader = csv.DictReader(csvfile)
        colum = [row[index] for row in reader]
        return np.array(colum)
def mainilearn(name,feaPse,TestPosName):
    for i in range(len(feaPse)):
        fea=feaPse[i]
        fileTrain = name+"model/"+ fea + ".csv"
        fileOriTest = name+"temp/"  + TestPosName[:-6]+fea + ".csv"
        outfileTest =  name+"temp/"  + TestPosName[:-6]+fea + "MRMD.csv"
        rows = colNum2(fileTrain)
        matrixArray = colReal2(fileOriTest, index=rows[0])
        for i in range(len(rows)):
            if i==0:
                continue
            else:
                matrix = colReal2(fileOriTest, index=rows[i])
                matrixArray = np.row_stack((matrixArray, matrix))
        realMatrix = matrixArray.T
        # np.savetxt(outfile, realMatrix, delimiter=",")
        pd_data = pd.DataFrame(realMatrix)
        pd_data.to_csv(outfileTest)
        print(fea+"end....")

def merge2Test(file0,file1,file2,file3,file4,file5,file6,file7,file8,file9,outfileArff):
    matrix0 = np.loadtxt(open(file0, "rb"), delimiter=",", skiprows=1)
    matrix0 = np.delete(matrix0, 0, axis=1)
    matrix1 = np.loadtxt(open(file1, "rb"), delimiter=",", skiprows=1)
    matrix1=np.delete(matrix1,0,axis=1)
    matrix2 = np.loadtxt(open(file2, "rb"), delimiter=",", skiprows=1)
    matrix2 = np.delete(matrix2, 0, axis=1)
    matrix3 = np.loadtxt(open(file3, "rb"), delimiter=",", skiprows=1)
    matrix3 = np.delete(matrix3, 0, axis=1)
    matrix4 = np.loadtxt(open(file4, "rb"), delimiter=",", skiprows=1)
    matrix4 = np.delete(matrix4, 0, axis=1)
    matrix5 = np.loadtxt(open(file5, "rb"), delimiter=",", skiprows=1)
    matrix5 = np.delete(matrix5, 0, axis=1)
    matrix6 = np.loadtxt(open(file6, "rb"), delimiter=",", skiprows=1)
    matrix6 = np.delete(matrix6, 0, axis=1)
    matrix7 = np.loadtxt(open(file7, "rb"), delimiter=",", skiprows=1)
    matrix7 = np.delete(matrix7, 0, axis=1)
    matrix8 = np.loadtxt(open(file8, "rb"), delimiter=",", skiprows=1)
    matrix8 = np.delete(matrix8, 0, axis=1)
    matrix9 = np.loadtxt(open(file9, "rb"), delimiter=",", skiprows=1)
    matrix9 = np.delete(matrix9, 0, axis=1)
    print(matrix0.shape)
    print(matrix1.shape)
    print(matrix2.shape)
    print(matrix3.shape)
    print(matrix4.shape)
    print(matrix5.shape)
    print(matrix6.shape)
    print(matrix7.shape)
    print(matrix8.shape)
    print(matrix9.shape)
    matrixTemp0 = np.hstack((matrix0, matrix1))
    matrixTemp1 = np.hstack((matrixTemp0, matrix2))
    matrixTemp2 = np.hstack((matrixTemp1, matrix3))
    matrixTemp3 = np.hstack((matrixTemp2, matrix4))
    matrixTemp4 = np.hstack((matrixTemp3, matrix5))
    matrixTemp5 = np.hstack((matrixTemp4, matrix6))
    matrixTemp6 = np.hstack((matrixTemp5, matrix7))
    matrixTemp7 = np.hstack((matrixTemp6, matrix8))
    matrix = np.hstack((matrixTemp7, matrix9))
    csv2arffTest(matrix,outfileArff)
    print("AllFea_end...")

def csv2arffTest(matrix,outfile):
    outfile1 = outfile[:-4] + ".csv"
    pd_data = pd.DataFrame(matrix)
    pd_data.to_csv(outfile1,header=None,columns=None)

def mergeSeq(name,TestPosName):
    pathTest = name + "temp/"
    file0 = "Kmer"
    file1 = "CKSNAP"
    file2 = "DNC"
    file3 = "Mismatch"
    file4 = "PC-PseDNC"
    file5 = "PC-PseTNC"
    file6 = "SC-PseDNC"
    file7 = "SC-PseTNC"
    file8 = "Fea"
    file9="NCP"
    outfileArffTest1 = pathTest + TestPosName[:-6]+"AllMRMD.csv"
    print(outfileArffTest1)
    merge2Test(pathTest + TestPosName[:-6]+ file0 + "MRMD.csv", pathTest + TestPosName[:-6]+ file1 + "MRMD.csv",
               pathTest +TestPosName[:-6] + file2 + "MRMD.csv",
               pathTest + TestPosName[:-6] + file3 + "MRMD.csv",
               pathTest + TestPosName[:-6]+ file4 + "MRMD.csv", pathTest + TestPosName[:-6] + file5 + "MRMD.csv",
               pathTest + TestPosName[:-6] + file6 + "MRMD.csv",
               pathTest + TestPosName[:-6] + file7 + "MRMD.csv", pathTest +TestPosName[:-6] + file8 + "MRMD.csv",
               pathTest + TestPosName[:-6] + file9 + ".csv",
               outfileArff=outfileArffTest1)

def Xg(modelSite,X_indipendent):


    from numpy import loadtxt
    import xgboost
    import pickle


    loaded_model = pickle.load(open(modelSite, "rb"))
    y_indipendent_pred=loaded_model.predict(X_indipendent)
    y_indipendent_prob = loaded_model.predict_proba(X_indipendent)
    indipendent_confidence = y_indipendent_prob[:, 1]
    return y_indipendent_pred,indipendent_confidence

def allFile2(path):
    import os
    files = os.listdir(path)
    f=[]
    for file in files:
        if ".dat" in file:
            f.append(path+file)
    return f

def getMatrix(fileTrain):
    X_train = np.loadtxt(open(fileTrain, "rb"), delimiter=",", skiprows=0)
    X_train = np.delete(X_train, 0, axis=1)
    return X_train

def seq():
    n = 41
    starting_point = 0
    samples = [sequence[i:i + n] for i in range(starting_point, len(sequence), n)]
    print(samples)

def seqfunction(sequence):
    n = 41
    starting_point = 0
    samples=[sequence[i:i+n] for i in range(starting_point, len(sequence), n)]
    return samples

def fasta(path,file):
    f = open(path+"Data/"+file)
    seq = {}
    for line in f:
        if line.startswith('>'):
            name = line.replace('>', '').split()[0]
            seq[name] = ''
        else:
            seq[name] += line.replace('\n', '').strip()
    f.close()
    return seq

def preproces(path,file):
    seq = fasta(path, file)
    values = list(seq.values())
    keys = list(seq.keys())
    fw = open(path +"Data/"+ file[:-6] + "Split.fasta", "w")
    seqName=[]
    for i in range(len(keys)):
        samples = seqfunction(values[i])
        for sample in samples:
            if len(sample) == 41:
                if "GGACA" or "GGACC" or "GGACU" or "GAACA" or "GAACC" or "GAACU" or "AGACA" or "AGACC" or "AGACU" or "AAACA" or "AAACC" or "AAACU" in sample:
                    fw.write(">"+keys[i] + "\n")
                    seqName.append(keys[i])
                    fw.write(sample + "\n")
    return seqName

def geneTxt2Csv(path,path2,trainfile):
    fr=open(path+trainfile,"r")
    line=fr.readlines()
    header=line[0]
    header=header.replace("\n","")
    header=header.replace("\"","")
    header = header.replace("\"", "")
    header=header.split(" ")
    header.insert(0,"class")
    matrix=[]
    for line in line[1:]:
        line=line.replace("FALSE","0")
        line=line.replace("TRUE", "1")
        line=line.split(" ")
        line=line[1:]
        line = list(map(float, line))
        matrix.append(line)
    trainDF=np.array(matrix)
    trainDF.to_csv(path2 + trainfile[:-7] + "Fea.csv", index=False, encoding='utf-8', header=False)
    print(trainfile+"train_end.......")

def GeneTxt2Csv(path,trainfile):
    # def geneTxt2Csv(path, path2, trainfile, trainpos, trainneg):
    fr = open(path + trainfile, "r")
    line = fr.readlines()
    header = line[0]
    header = header.replace("\n", "")
    header = header.replace("\"", "")
    header = header.replace("\"", "")
    header = header.split(" ")
    header.insert(0, "class")
    matrix = []
    for line in line[1:]:
        line = line.replace("FALSE", "0")
        line = line.replace("TRUE", "1")
        line = line.split(" ")
        line = line[1:]
        line = list(map(float, line))
        matrix.append(line)
    train = np.array(matrix)

    classCol = np.hstack((np.ones(1), np.zeros(train.shape[0]-1)))
    print(classCol.shape)
    train = np.column_stack((classCol, train))
    train = np.vstack((header, train))
    # print(train)
    trainDF = pd.DataFrame(train)
    trainDF.to_csv(path + trainfile[:-4] + ".csv", index=False, encoding='utf-8', header=False)
    print(trainfile + "end.......")
import os
def main(TestPosName):

    name = os.getcwd()+"/"
    seqName=preproces(name,TestPosName)
    TestPosName=TestPosName[:-6] + "Split.fasta"
    ############################ Data
    ####################### featurextract
    main1(name,TestPosName)
    main2(name,TestPosName)
    main3(name,TestPosName)
    main4(name, TestPosName)
    orderPse(name,TestPosName)
    txt2csv(name,TestPosName)
    os.system("Rscript Fea.R")
    GeneTxt2Csv(name+"temp/", TestPosName[:-6]+"Fea.txt")
    ##################################feaselection
    feas = ["Kmer", "DNC", "CKSNAP", "Mismatch", "PC-PseDNC", "PC-PseTNC", "SC-PseDNC", "SC-PseTNC","Fea"]
    mainilearn(name, feas,TestPosName)
    mergeSeq(name, TestPosName)
    X_test = getMatrix(name+"temp/" +TestPosName[:-6]+ "AllMRMD.csv")
    path = name+"model/MatureModel/"
    model = allFile2(path)
    labelList = []
    confList = []
    all = []
    for model in model:
        print(model)
        label, conf = Xg(model, X_test)
        labelList.append(label)
        confList.append(conf)
        all.append(label)
        all.append(conf)
    allArray = np.array(all)
    seqName = np.array(seqName)
    seqAll = np.vstack((seqName, allArray))
    pre=[]
    for i in range(len(confList[0])):
        if (confList[0][i]+confList[1][i]+confList[2][i]+confList[3][i]+confList[4][i]+confList[5][i])>0.5:
            pre.append("M6A")
        else:
            pre.append("non_M6A")
    pre = []
    for i in range(len(confList[0])):
        if (confList[0][i] + confList[1][i] + confList[2][i] + confList[3][i] + confList[4][i] + confList[5][i]) / 6 > 0.5:
            pre.append("M6A")
        else:
            pre.append("non_M6A")
    pre2 = np.array(pre)
    pre = np.vstack((seqAll, pre2))
    conf = np.mean(np.array(confList), axis=0)
    a = np.vstack((pre, conf)).T
    header = ["seqName", "A549Label", "A549Conf", "CD8TLabel", "CD8TConf", "HEK_abacmLabel", "HEK_abacmConf",
              "HEK_sysyLabel", "HEK_sysyConf", "HeLaLabel", "HeLaConf", "MOLMLabel", "MOLMConf", "voteLabel",
              "voteConf"]
    c = np.vstack((header, a))
    b = pd.DataFrame(c)
    b.to_csv(name + "Result/" + TestPosName[:-6] + 'Mature.csv', index=False, encoding='utf-8', header=None)
    # os.system("rm -r /home/lijing/HSM6AP/temp/*")
import argparse
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-infa", "--in_fa", action="store", dest='in_fa', required=True,
                        help="input fa")
    args = parser.parse_args()
    in_fa = args.in_fa
    main(in_fa)


