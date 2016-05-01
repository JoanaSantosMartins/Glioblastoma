# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:48:47 2016

@author: Joana
"""

class netMHCpan():
    
#Objetivo do Script:
#Processar os outputs do netMHCpan e netMHCIIpan
    
    
    def __init__(self, filename = []):
        self.filename = filename
    
    def edit(self, feature):
        #Função que permite editar/processar texto proveniente de ficheiros.
        
        f = str(feature)
        f = f.replace("[", "")
        f = f.replace("]", "")
        f = f.replace("'", "")
        f = f.replace(" ", "")
        return f

    def remove_head(self):
        #Função que remove os cabeçalhos.
    
        file = open(self.filename, "r")
        lines = file.read().splitlines(True)
        new_file = open("u1.txt", "w")
        new_file.writelines(lines[50:])
        new_file.close() 
        file.close()

    def remove_lines(self):
        #Função que remove linhas com informação irrelevante para a análise.        
        
        inFile = open(self.filename)
        outFile = open("result.txt", "w")
        save = []
        keepSave = True
        for line in inFile:
            save.append(line)
            if line.startswith("-"):
                if keepSave:
                    outFile.write("".join(save))
                keepSave = False
                save = []
            elif line.startswith("  p"):
                if keepSave:
                    outFile.write("".join(save))
                keepSave = False
                save = []
            elif line.startswith("  "):
                keepSave = True
        inFile.close()
        outFile.close()

        
    def human_HLA(self):
        #Função que permite obter uma lista com os HLA numa única coluna.        
        
        inFile = open(self.filename)
        outFile = open("HLA.txt", "w")
        save = []
        for line in inFile:
            save.append(line)
        save = [i.replace("\n","") for i in save]
        outFile.write(str(self.edit(save)))
        inFile.close()
        outFile.close()
        
    def read_output_I(self):
        #Função que lê todos os outputs do netMHCpan e guarda num ficheiro todas as linhas correspondentes a SB e WB.
    
        import glob   
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/out/*.txt'   
        files=glob.glob(path) 
        for file in files:  
            f=open(file, 'r') 
            lines = f.readlines()
            for line in lines:
                line = line.strip()           
                column = line.split()
                if column[8:9] == ['WB'] or column[8:9] == ['SB']:
                    outFile = open("14_SB_WB.txt", "a")
                    outFile.write(str(line) + "\n")
                    outFile.close()
            f.close()
            
    def read_output_II(self):
        #Função que lê todos os outputs do netMHCIIpan e guarda num ficheiro todas as linhas correspondentes a SB e WB.
    
        import glob   
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/out/*.txt'   
        files=glob.glob(path) 
        for file in files:  
            f=open(file, 'r') 
            lines = f.readlines()
            for line in lines:
                line = line.strip()           
                column = line.split()

                if column[11:12] == ['<=WB'] or column[11:12] == ['<=SB']:
                    outFile = open("9_SB_WB.txt", "a")
                    outFile.write(str(line) + "\n")
                    outFile.close()

            f.close()
        
if __name__=="__main__":
    
    def test1(): 
        n = netMHCpan("result.txt")
        n.remove_head()
        n.remove_lines()
    
    def test2():
        n = netMHCpan("human_HLA.txt")
        n.human_HLA()
        
    def test3():
        n = netMHCpan()
        n.read_output_II()
    
    
    test3()
        
    
    
