# -*- coding: utf-8 -*-
"""
Created on Sun May 22 17:24:53 2016

@author: Joana
"""

import glob 
import re

class antigenicity():
#Objetivo do Script:
#Avaliar antigenicidade da proteína
        
    def __init__(self, n_mut, length, pos):
        self.n_mut = n_mut #número de mutações
        self.length = length #comprimento da sequência
        self.pos = pos #posição da mutação
        self.mut_SB = self.create_zeros(self.n_mut, self.length) #analisadas diferentes mutações
        self.mut_WB = self.create_zeros(self.n_mut, self.length)
        self.wt_SB = [0] * self.length #analisada toda a sequência
        self.wt_WB = [0] * self.length
    
    def create_zeros(self, lin, col):
        #Função que cria matriz de zeros (número de mutações x tamanho da sequência).
        
        mat = []
        for i in range(0, lin):
            mat.append([0]*col)
        return mat

    def transposta(self, mat):
        #Função que retorna a transposta de uma matriz.
    
        t = self.create_zeros(self.length, self.n_mut)
        for i in range(0,self.n_mut):
            for j in range(0,self.length):
                t[j][i] = mat[i][j]
        return t    
    
    def edit(self, feature):
        #Função que permite editar/processar texto proveniente de ficheiros.
        
        f = str(feature)
        f = f.replace("[", "")
        f = f.replace("]", "")
        f = f.replace("*", "")
        f = f.replace("'", "")
        f = f.replace(", ", " ")
        return f
    
    def count_mut(self, path_pep):  
        #Função que faz a contagem dos WB e SB para cada aminoácido na sequência com a mutação.        
        
        files = glob.glob(path_pep) 
        for file in files: 
            f = open(file, "r") 
            lines = f.readlines()
            for line in lines:
                line = line.strip()           
                column = line.split()
                len_peptide = len(self.edit(column[2:3]))
                pat = 'MUT[0-9]'
                m = re.search(pat, line)    
                
                if m and (column[8:9] == ['SB'] or column[11:12] == ['<=SB']):
                    for i in self.edit(column[3:4]):
                        
                        if i.isdigit():   
                            i = int(i) - 1
                            ini = self.pos[i] - len_peptide + int(self.edit(column[0:1])) 
                            fin = self.pos[i] + int(self.edit(column[0:1]))                            
                            for j in range(ini, fin):
                                self.mut_SB[i][j] += 1
 
                elif column[3:4] == ['MUT'] and (column[8:9] == ['SB'] or column[11:12] == ['<=SB']):
                    ini = self.pos[0] - len_peptide + int(self.edit(column[0:1])) 
                    fin = self.pos[0] + int(self.edit(column[0:1]))                            
                    for j in range(ini, fin):
                        self.mut_SB[0][j] += 1
               
                if m and (column[8:9] == ['WB'] or column[11:12] == ['<=WB']):
                    for i in self.edit(column[3:4]):
                        if i.isdigit():   
                            i = int(i) - 1
                            ini = self.pos[i] - len_peptide + int(self.edit(column[0:1])) 
                            fin = self.pos[i] + int(self.edit(column[0:1]))                            
                            for j in range(ini, fin):
                                self.mut_WB[i][j] += 1
                elif column[3:4] == ['MUT'] and (column[8:9] == ['WB'] or column[11:12] == ['<=WB']):
                    ini = self.pos[0] - len_peptide + int(self.edit(column[0:1])) 
                    fin = self.pos[0] + int(self.edit(column[0:1]))                            
                    for j in range(ini, fin):
                        self.mut_WB[0][j] += 1
            f.close()
    
    def count_wt(self, path_seq):
        #Função que faz a contagem dos WB e SB para cada aminoácido na sequência wild-type. 
        
        files = glob.glob(path_seq)        
        for file in files:  
            f = open(file, "r") 
            lines = f.readlines()
            for line in lines:
                line = line.strip()           
                column = line.split()
                len_peptide = len(self.edit(column[2:3]))
                
                if column[8:9] == ['SB'] or column[11:12] == ['<=SB']:
                    ini = int(self.edit(column[0:1])) 
                    fin = len_peptide + int(self.edit(column[0:1]))
                    for i in range(ini, fin):
                        self.wt_SB[i] += 1
                     
                if column[8:9] == ['WB'] or column[11:12] == ['<=WB']:
                    ini = int(self.edit(column[0:1])) 
                    fin = len_peptide + int(self.edit(column[0:1])) 
                    for i in range(ini, fin):                      
                        self.wt_WB[i] += 1
                
            f.close()

    def bind(self):
        #Função que substituiu os zeros da sequência com mutação pelo respetivo número de SB  e WB.

        for i in range(0, self.n_mut):
            for j in range(0, self.length):                
                if self.mut_SB[i][j] == 0:
                    self.mut_SB[i][j] = self.wt_SB[j]
                if self.mut_WB[i][j] == 0:
                    self.mut_WB[i][j] = self.wt_WB[j]

    def save_file(self, file_name):
        #Função que guarda dos valores de cada lista num ficheiro, dispondo os valores em colunas.
        
        imuno = open(file_name, "w") 
        for i in range(0, self.n_mut):
            j = i + 1
            imuno.write("WB_Mut" + str(j) + " ")
        imuno.write("WB_Wt ")
        for i in range(0, self.n_mut):
            j = i + 1
            imuno.write("SB_Mut" + str(j) + " ")
        imuno.write("SB_Wt " + "\n")
        wb = self.transposta(self.mut_WB)
        sb = self.transposta(self.mut_SB)
        for i in range(0, self.length):
            imuno.write(str(self.edit(wb[i])) + " " + str(self.wt_WB[i]) + " " + str(self.edit(sb[i])) + " " + str(self.wt_SB[i]) + "\n") 
        imuno.close()



if __name__=="__main__":
    
    def test1(): #IDH1
        path_pep = "D:/PI/Bioinformatics/2semestre/Projeto/scripts/IDH1/SB_WB_II/pep/*.txt"
        path_seq = "D:/PI/Bioinformatics/2semestre/Projeto/scripts/IDH1/SB_WB_II/seq/*.txt"
        i = antigenicity(1, 414, [132])
        i.count_mut(path_pep)
        i.count_wt(path_seq)
        i.bind()
        i.save_file("imuno_idh1_II.txt")
    
    def test2(): #IDH2
        path_pep = "D:/PI/Bioinformatics/2semestre/Projeto/scripts/IDH2/SB_WB_II/pep/*.txt"
        path_seq = "D:/PI/Bioinformatics/2semestre/Projeto/scripts/IDH2/SB_WB_II/seq/*.txt"
        i = antigenicity(3, 452, [172,172,172])
        i.count_mut(path_pep)
        i.count_wt(path_seq)
        i.bind()
        i.save_file("imuno_idh2_II.txt")
    
    def test3(): #TP53
        path_pep = "D:/PI/Bioinformatics/2semestre/Projeto/scripts/TP53/SB_WB_II/pep/*.txt"
        path_seq = "D:/PI/Bioinformatics/2semestre/Projeto/scripts/TP53/SB_WB_II/seq/*.txt"
        i = antigenicity(5, 393, [143,175,248,273,281])
        i.count_mut(path_pep)
        i.count_wt(path_seq)
        i.bind()
        i.save_file("imuno_tp53_II.txt")
    
    def test4(): #Histone H3-3
        path_pep = "D:/PI/Bioinformatics/2semestre/Projeto/scripts/H3/SB_WB_II/pep/*.txt"
        path_seq = "D:/PI/Bioinformatics/2semestre/Projeto/scripts/H3/SB_WB_II/seq/*.txt"
        i = antigenicity(3, 136, [28,35,35])
        i.count_mut(path_pep)
        i.count_wt(path_seq)
        i.bind()
        i.save_file("imuno_h3_II.txt")


    test1()
    test2()
    test3()
    test4()