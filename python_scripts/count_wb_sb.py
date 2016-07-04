# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 00:46:25 2016

@author: Joana
"""

import glob  
import re

class count_wb_sb():
    
#Objetivo do Script:
#Contagem de SB e WB para wildtype e mutante.
    
    def __init__(self, file, n_mut):
        self.file = file
        self.n_mut = n_mut #número de mutações
        self.alleles = [] #lista com HLA class I ou class II
        self.count = {} #dicionário com HLA {nome HLA: SB WT, WB WT, SB MUT, WB MUT}
        self.mut_SB = self.human_HLA()
        self.mut_WB = self.human_HLA()
        self.wt_SB = self.human_HLA()
        self.wt_WB = self.human_HLA()
    
    def create_zeros(self, lin, col):
        #Função que cria matriz de zeros.
        
        mat = []
        for i in range(0, lin):
            mat.append([]*col)
        return mat
    
    def edit(self, feature):
        #Função que permite editar/processar texto proveniente de ficheiros.
        
        f = str(feature)
        f = f.replace("[", "")
        f = f.replace("]", "")
        f = f.replace("*", "")
        f = f.replace("'", "")
        f = f.replace(", ", " ")
        return f    
    
    def human_HLA(self):
        #Função que lê o ficheiro com os HLA e guarda num dicionário.        
        
        res = {}
        f = open(self.file)
        for line in f:
            self.alleles.append(line)
        self.alleles = [i.replace("\n","") for i in self.alleles]
        res = dict([(x,[0]*self.n_mut) for x in self.alleles])
        f.close()  
        return res

    def SB_WB(self, path):
        #Função que faz a contagem de WB e SB para cada HLA.        
        
        files = glob.glob(path) 
        for file in files:  
            f = open(file, "r") 
            lines = f.readlines()
            for line in lines:
                line = line.strip()           
                column = line.split()
                pat1 = 'MUT[0-9]'
                pat2 = 'WT[0-9]'
                m1 = re.search(pat1, line)
                m2 = re.search(pat2, line)
                
                if m1:
                    for i in self.edit(column[3:4]):
                        if i.isdigit():
                            i = int(i) - 1
                            if column[8:9] == ['SB'] or column[11:12] == ['<=SB']: 
                                self.mut_SB[str(self.edit(column[1:2]))][i] += 1
                            elif column[8:9] == ['WB'] or column[11:12] == ['<=WB']: 
                                self.mut_WB[str(self.edit(column[1:2]))][i] += 1
                elif m2:
                    for i in self.edit(column[3:4]):
                        if i.isdigit():
                            i = int(i) - 1
                            if column[8:9] == ['SB'] or column[11:12] == ['<=SB']: 
                                self.wt_SB[str(self.edit(column[1:2]))][i] += 1
                            elif column[8:9] == ['WB'] or column[11:12] == ['<=WB']: 
                                self.wt_WB[str(self.edit(column[1:2]))][i] += 1
                                
                else:
                    if (column[8:9] == ['SB'] and column[3:4] == ['WT']) or (column[11:12] == ['<=SB'] and column[3:4] == ['WT']): 
                        self.wt_SB[str(self.edit(column[1:2]))][0] += 1
                    elif (column[8:9] == ['WB'] and column[3:4] == ['WT']) or (column[11:12] == ['<=WB'] and column[3:4] == ['WT']): 
                        self.wt_WB[str(self.edit(column[1:2]))][0] += 1
                    elif (column[8:9] == ['SB'] and column[3:4] == ['MUT']) or (column[11:12] == ['<=SB'] and column[3:4] == ['MUT']): 
                        self.mut_SB[str(self.edit(column[1:2]))][0] += 1
                    elif (column[8:9] == ['WB'] and column[3:4] == ['MUT']) or (column[11:12] == ['<=WB'] and column[3:4] == ['MUT']): 
                        self.mut_WB[str(self.edit(column[1:2]))][0] += 1
            f.close()  
                    
    def total_WB_SB(self,file):
        #Função que faz a contagem de todos os SB e WB para Wt e Mut.
        
        WB = self.create_zeros(self.n_mut, 2) #wt, mut
        SB = self.create_zeros(self.n_mut, 2) #wt, mut 
        for key1, value1 in self.mut_SB.items():
            for key2, value2 in self.wt_SB.items():
                if key1 == key2:
                    for i in range(0, self.n_mut):
                        SB[i][0] += value2[i]
                        SB[i][1] += value1[i]
        for key1, value1 in self.mut_WB.items():
            for key2, value2 in self.wt_WB.items():
                if key1 == key2:
                    for i in range(0, self.n_mut):
                        WB[i][0] += value2[i]
                        WB[i][1] += value1[i]      
        f = open(file, "w") 
        for i in range(0, self.n_mut):
            j = i + 1
            f.write("Wt" + str(j) + " " + "Mut" + str(j) + " ")
        f.write("\n" + str(self.edit(WB))  + "\n" + str(self.edit(SB)))
        f.close()
    
    def best_HLA(self):
        #Função que faz o ranking dos HLA com maior impacto.
        
        SB = self.create_zeros(self.n_mut, 1)
        WB = self.create_zeros(self.n_mut, 1)
        
        for key1, value1 in self.mut_SB.items():
            for key2, value2 in self.wt_SB.items():
                if key1 == key2:
                    for i in range(0, self.n_mut):  
                            if value1[i] > value2[i]:
                                diff = value1[i] - value2[i]
                                SB[i].append((key1, value1[i], diff))
        for i in SB:
            for j in SB:
                for tup in j:
                    j.sort(key=lambda tup:(tup[2], tup[1]), reverse = True)
        
                
        for key1, value1 in self.mut_WB.items():
            for key2, value2 in self.wt_WB.items():
                if key1 == key2:
                    for i in range(0, self.n_mut):
                        if value1[i] > value2[i]:
                            diff = value1[i] - value2[i]
                            WB[i].append((key1, value1[i], diff)) 
        
        for i in WB:
            for j in WB:
                for tup in j:
                    j.sort(key=lambda tup:(tup[2], tup[1]), reverse = True)
        
        return SB
        



if __name__=="__main__":
    
    #IDH1
    def test1(): #netMHCpan 
        file = "alleles_HLA_I.txt"
        c = count_wb_sb(file, 1)
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/IDH1/SB_WB_I/pep/*.txt'
        c.SB_WB(path)
        #c.total_WB_SB("sb_wb_idh1_I.txt")
        c.best_HLA()
        
    def test2(): #netMHCIIpan 
        file = "alleles_HLA_II.txt"
        c = count_wb_sb(file, 1)
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/IDH1/SB_WB_II/pep/*.txt'
        c.SB_WB(path)
        #c.total_WB_SB("sb_wb_idh1_II.txt")
        c.best_HLA()
    
    #IDH2
    def test3(): #netMHCpan 
        file = "alleles_HLA_I.txt"
        c = count_wb_sb(file, 3)
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/IDH2/SB_WB_I/pep/*.txt'
        c.SB_WB(path)
        #c.total_WB_SB("sb_wb_idh2_I.txt")
        c.best_HLA()
        
    def test4(): #netMHCIIpan 
        file = "alleles_HLA_II.txt"
        c = count_wb_sb(file, 3)
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/IDH2/SB_WB_II/pep/*.txt'
        c.SB_WB(path)
        #c.total_WB_SB("sb_wb_idh2_II.txt")
        c.best_HLA()
    
    #TP53
    def test5(): #netMHCpan 
        file = "alleles_HLA_I.txt"
        c = count_wb_sb(file, 5)
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/TP53/SB_WB_I/pep/*.txt'
        c.SB_WB(path)
        #c.total_WB_SB("sb_wb_tp53_I.txt")
        c.best_HLA()
        
    def test6(): #netMHCIIpan 
        file = "alleles_HLA_II.txt"
        c = count_wb_sb(file, 5)
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/TP53/SB_WB_II/pep/*.txt'
        c.SB_WB(path)
        #c.total_WB_SB("sb_wb_tp53_II.txt")
        c.best_HLA()
    
    #H3
    def test7(): #netMHCpan 
        file = "alleles_HLA_I.txt"
        c = count_wb_sb(file, 3)
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/H3/SB_WB_I/pep/*.txt'
        c.SB_WB(path)
        #c.total_WB_SB("sb_wb_h3_I.txt")
        c.best_HLA()
        
    def test8(): #netMHCIIpan 
        file = "alleles_HLA_II.txt"
        c = count_wb_sb(file, 3)
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/H3/SB_WB_II/pep/*.txt'
        c.SB_WB(path)
        #c.total_WB_SB("sb_wb_h3_II.txt")
        c.best_HLA()
    
    test8() 



