# -*- coding: utf-8 -*-
"""
Created on Wed Mar 23 15:48:47 2016

@author: Joana
"""

import glob

class netMHCpan():
    
#Objetivo do Script:
#Processar os outputs do netMHCpan e netMHCIIpan
    
    
    def __init__(self, path, length, netMHC):
        self.path = path #diretoria
        self.length = length #comprimento do peptido
        self.netMHC = netMHC #netMHCpan = 1 e netMHCIIpan = 2
    
    def edit(self, feature):
        #Função que permite editar/processar texto proveniente de ficheiros.
        
        f = str(feature)
        f = f.replace("[", "")
        f = f.replace("]", "")
        f = f.replace("'", "")
        f = f.replace(" ", "")
        return f
        
    def remove_lines(self):
        #Função que remove linhas com informação irrelevante para a análise.         
        
        files = glob.glob(self.path) 
        for file in files:  
            f = open(file, "r") 
            lines = f.readlines()
            f.close()   
            f2 = open(file, "w")
            for line in lines:
                line = line.strip()           
                column = line.split()
                pos = self.edit(column[0:1])
                if pos.isdigit():
                    f2.write(line + "\n")
            f2.close()
        
    def remove_error(self):
        #Função que remove erros nos péptidos.
        
        files = glob.glob(self.path) 
        for file in files:  
            f = open(file, "r") 
            lines = f.readlines()
            f.close()   
            f2 = open(file, "w")
            for line in lines:
                line = line.strip()           
                column = line.split()
                
                if self.netMHC == 1: #apresenta péptidos com aminoácido X.
                    pep = self.edit(column[2:3])
                    error = False
                    for a in pep:
                        if a == "X":
                            error = True
                    if error == False:
                            f2.write(line + "\n")
                            
                if self.netMHC == 2: #apresenta péptidos com tamanho inferior.
                    if len(self.edit(column[2:3])) == self.length:
                        f2.write(line + "\n")
            f2.close()
        
    def find_SB_WB(self):
        #Função que lê todos os outputs do netMHCpan e guarda num ficheiro todas as linhas correspondentes a SB e WB.
    
        files = glob.glob(self.path) 
        for file in files:  
            f = open(file, "r") 
            lines = f.readlines()
            for line in lines:
                line = line.strip()           
                column = line.split()
                
                if self.netMHC == 1:
                    if column[8:9] == ['WB'] or column[8:9] == ['SB']:
                        outFile = open(str(self.length) + "_SB_WB.txt", "a")
                        outFile.write(str(line) + "\n")
                        outFile.close()
                    
                if self.netMHC == 2:
                    if column[11:12] == ['<=WB'] or column[11:12] == ['<=SB']:
                        outFile = open(str(self.length) + "_SB_WB.txt", "a")
                        outFile.write(str(line) + "\n")
                        outFile.close()
            f.close()
            
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
        
if __name__=="__main__":
    
    def test1(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCpan-2.8/H3/seq/out9/*.txt"
        n = netMHCpan(path, 9, 1)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test2(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCpan-2.8/H3/seq/out10/*.txt"
        n = netMHCpan(path, 10, 1)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test3(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCpan-2.8/H3/seq/out11/*.txt"
        n = netMHCpan(path, 11, 1)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test4(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCpan-2.8/H3/seq/out12/*.txt"
        n = netMHCpan(path, 12, 1)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test5(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCpan-2.8/H3/seq/out13/*.txt"
        n = netMHCpan(path, 13, 1)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test6(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCpan-2.8/H3/seq/out14/*.txt"
        n = netMHCpan(path, 14, 1)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test7(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCpan-2.8/H3/seq/out8/*.txt"
        n = netMHCpan(path, 8, 1)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test8(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCIIpan-3.1/seq/out16/*.txt"
        n = netMHCpan(path, 16, 2)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test9(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCIIpan-3.1/seq/out17/*.txt"
        n = netMHCpan(path, 17, 2)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test10(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCIIpan-3.1/seq/out18/*.txt"
        n = netMHCpan(path, 18, 2)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test11(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCIIpan-3.1/seq/out19/*.txt"
        n = netMHCpan(path, 19, 2)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
    
    def test12(): 
        path = "D:/PI/Bioinformatics/2semestre/Projeto/Glioblastoma-master_v2/netMHCIIpan-3.1/seq/out20/*.txt"
        n = netMHCpan(path, 20, 2)
        n.remove_lines()
        n.remove_error()
        n.find_SB_WB()
             
    def test15():
        n = netMHCpan("human_HLA.txt")
        n.human_HLA()

      
    test1()
    test2()    
    test3()   
    test4() 
    test5()
    
    test6()
    test7()

    
    
    