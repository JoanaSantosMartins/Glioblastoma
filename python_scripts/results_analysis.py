# -*- coding: utf-8 -*-
"""
Created on Sun Apr 24 00:46:25 2016

@author: Joana
"""
from openpyxl import Workbook
from collections import OrderedDict

class results_analysis():
    
#Objetivo do Script:
#Contagem de SB e WB para wildtype e mutante, guardar resultados num ficheiro excel.
    
    def __init__(self):
        self.alleles_I = [] #lista com HLA class I
        self.alleles_II = [] #lista com HLA class II
        self.count_I = {} #dicionário com HLA class I {nome HLA: SB WT, SB MUT, WB WT, WB MUT}
        self.count_II = {} #dicionário com HLA class II {nome HLA: SB WT, SB MUT, WB WT, WB MUT}
    
    def human_HLA(self):
        #Função que lê o ficheiro com os HLA e guarda num dicionário.        
        
        file_I = open("alleles_HLA_I.txt")
        file_II = open("alleles_HLA_II.txt")
        for line in file_I:
            self.alleles_I.append(line)
        self.alleles_I = [i.replace("\n","") for i in self.alleles_I]
        self.count_I = dict([(x,[int(0),int(0),int(0),int(0)]) for x in self.alleles_I])

        for line in file_II:
            self.alleles_II.append(line)
        self.alleles_II = [i.replace("\n","") for i in self.alleles_II]
        self.count_II = dict([(x,[int(0),int(0),int(0),int(0)]) for x in self.alleles_II])
        file_I.close()  
        file_II.close() 
    
    def main_xls(self, ws1, ws2):
        #Função que cria um ficheiro excel com a contagem do WB e SB.        
        
        ws1['A1'] = "HLA class I"
        ws1['B2'] = "Number of peptides" 
        ws1['A6'] = "Alleles"
        ws1['B6'] = "SB Wild-type"
        ws1['C6'] = "SB Mutant"
        ws1['E6'] = "WB Wild-type"
        ws1['F6'] = "WB Mutant"
        ws1['H1'] = "HLA class II"
        ws1['I2'] = "Number of peptides" 
        ws1['H6'] = "Alleles"
        ws1['I6'] = "SB Wild-type"
        ws1['J6'] = "SB Mutant"
        ws1['L6'] = "WB Wild-type"
        ws1['M6'] = "WB Mutant"

        ws2['A1'] = "HLA class I"
        ws2['B2'] = "Number of peptides" 
        ws2['A6'] = "Alleles"
        ws2['B6'] = "SB Wild-type"
        ws2['C6'] = "SB Mutant"
        ws2['E6'] = "WB Wild-type"
        ws2['F6'] = "WB Mutant"
        ws2['H1'] = "HLA class II"
        ws2['I2'] = "Number of peptides" 
        ws2['H6'] = "Alleles"
        ws2['I6'] = "SB Wild-type"
        ws2['J6'] = "SB Mutant"
        ws2['L6'] = "WB Wild-type"
        ws2['M6'] = "WB Mutant"
        
        I = OrderedDict(sorted(self.count_I.items(), key=lambda t: t[0]))
        II = OrderedDict(sorted(self.count_II.items(), key=lambda t: t[0]))
        #print(II)
        
        k = 6 
        for key, value in I.items():
            k += 1
            h = str(k)
            ws1['A' + h] = key
            ws1['B' + h] = value[0]
            ws1['C' + h] = value[2]
            ws1['E' + h] = value[1]
            ws1['F' + h] = value[3]
            
        k = 6        
        for key, value in I.items():
            if (value[0] <= value[2] and value[1] < value[3]):
                k += 1
                h = str(k)
                ws2['A' + h] = key
                ws2['B' + h] = value[0]
                ws2['C' + h] = value[2]
                ws2['E' + h] = value[1]
                ws2['F' + h] = value[3]
                
            if (value[0] <= value[2] and value[1] >= value[3]):
                k += 1
                h = str(k)
                ws2['A' + h] = key
                ws2['B' + h] = value[0]
                ws2['C' + h] = value[2]
            
            if (value[0] >= value[2] and value[1] < value[3]):
                k += 1
                h = str(k)
                ws2['A' + h] = key
                ws2['E' + h] = value[1]
                ws2['F' + h] = value[3]
            
            
        k = 6
        for key, value in II.items():
            k += 1
            h = str(k)
        
            ws1['H' + h] = key
            ws1['I' + h] = value[0]
            ws1['J' + h] = value[2]
            ws1['L' + h] = value[1]
            ws1['M' + h] = value[3]
        
        k = 6        
        for key, value in II.items():
            if (value[0] < value[2] and value[1] < value[3]):
                k += 1
                h = str(k)
                ws2['H' + h] = key
                ws2['I' + h] = value[0]
                ws2['J' + h] = value[2]
                ws2['L' + h] = value[1]
                ws2['M' + h] = value[3]
                
            if (value[0] < value[2] and value[1] >= value[3]):
                k += 1
                h = str(k)
                ws2['H' + h] = key
                ws2['I' + h] = value[0]
                ws2['J' + h] = value[2]

            if (value[0] >= value[2] and value[1] < value[3]):
                k += 1
                h = str(k)
                ws2['H' + h] = key
                ws2['L' + h] = value[1]
                ws2['M' + h] = value[3]
        

    
        return k
    
    def edit(self, feature):
        #Função que permite editar/processar texto proveniente de ficheiros.
        
        f = str(feature)
        f = f.replace("[", "")
        f = f.replace("]", "")
        f = f.replace("*", "")
        f = f.replace("'", "")
        f = f.replace(", ", "\n")
        return f
   
    def count_SB_WB_I(self):
        #Função que faz a contagem de WB e SB dos HLA de class I.        
        
        import glob   
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/SB_WB_I/*.txt'   
        files=glob.glob(path) 
        for file in files:  
            f=open(file, 'r') 
            lines = f.readlines()
            for line in lines:
                line = line.strip()           
                column = line.split()
                if str(self.edit(column[1:2])) in self.count_I.keys():
                    if column[8:9] == ['SB'] and column[3:4] == ['WT']: 
                        self.count_I[str(self.edit(column[1:2]))][0] += 1
                    elif column[8:9] == ['WB'] and column[3:4] == ['WT']: 
                        self.count_I[str(self.edit(column[1:2]))][1] += 1
                    elif column[8:9] == ['SB'] and column[3:4] == ['MUT']: 
                        self.count_I[str(self.edit(column[1:2]))][2] += 1
                    elif column[8:9] == ['WB'] and column[3:4] == ['MUT']: 
                        self.count_I[str(self.edit(column[1:2]))][3] += 1 
            f.close()     
    
    def count_SB_WB_II(self):
        #Função que faz a contagem de WB e SB dos HLA de class II.          
        
        import glob   
        path = 'D:/PI/Bioinformatics/2semestre/Projeto/scripts/SB_WB_II/*.txt'   
        files=glob.glob(path) 
        for file in files:  
            f=open(file, 'r') 
            lines = f.readlines()
            for line in lines:
                line = line.strip()           
                column = line.split()
                if str(self.edit(column[1:2])) in self.count_II.keys():
                    if column[11:12] == ['<=SB'] and column[3:4] == ['WT']:
                        self.count_II[str(self.edit(column[1:2]))][0] += 1
                    elif column[11:12] == ['<=WB'] and column[3:4] == ['WT']:  
                        self.count_II[str(self.edit(column[1:2]))][1] += 1
                    elif column[11:12] == ['<=SB'] and column[3:4] == ['MUT']: 
                        self.count_II[str(self.edit(column[1:2]))][2] += 1
                    elif column[11:12] == ['<=WB'] and column[3:4] == ['MUT']: 
                        self.count_II[str(self.edit(column[1:2]))][3] += 1

            f.close()  


if __name__=="__main__":
    
    def test1(): 
        ra = results_analysis()
        ra.human_HLA()
        ra.count_SB_WB_I()
        ra.count_SB_WB_II()
                
        
        wb = Workbook()
        filename = "IDH1_analysis.xlsx"
        ws1 = wb.active
        ws1.title = "IDH1_analysis"
        ws2 = wb.create_sheet(title = "MUT_SB_WB")
        ra.main_xls(ws1, ws2)
        wb.save(filename)
        


    test1()
        



