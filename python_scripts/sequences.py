# -*- coding: utf-8 -*-
"""
Created on Mon Mar 21 11:13:22 2016

@author: Joana
"""

from Bio import Entrez
from Bio import SeqIO
Entrez.email = "joana_smar@hotmail.com"


class sequence():

#Objetivo do Script:
#Obter a sequência da proteína com recurso ao NCBI
#Guardar a sequência em formato FASTA
#Criar péptidos para serem analisados pelo netMHCpan
    
    def __init__(self, id_protein, file_protein = None):
        self.id_protein = id_protein #id da proteína
        self.file_protein = file_protein
        self.sequence = "" #sequência da proteína
        self.mut_sequence = [] #sequências com mutações 
        self.definition = "" #informação sobre o nome da proteína e respetivo organismo associado

#Procurar a sequência por ID da proteína:

    def protein_ncbi(self):
        #Função que acede ao NCBI a partir do id de uma proteína, retirando a respetiva sequência.

        handle = Entrez.efetch(db = "protein", rettype = "gb", retmode = "xml", id = self.id_protein)
        results = Entrez.read(handle)
        self.sequence = results[0]["GBSeq_sequence"].upper()      
        self.definition = results[0]["GBSeq_definition"]
        return self.sequence, self.definition

#Obter a sequência a partir de um ficheiro fasta:

    def protein_file(self):
        #Função que obtém a sequência de uma proteína a partir de um ficheiro em formato fasta.        
        
        handle = open(self.file_protein, "r")
        for record in SeqIO.parse(handle, "fasta"):
            self.sequence = record.seq 
        handle.close()

#Guardar a sequência num ficheiro:
    
    def protein_fasta(self, protein_name):
        #Função que cria um ficheiro em formato fasta com a sequência da proteína obtida anteriormente.       
        
        file = open(protein_name + ".fasta", "w")
        file.write(">gi|" + self.id_protein + "| " + self.definition +  "\n" + self.sequence)
        file.close() 
        
#Induzir as mutações:

    def induce_mutation(self, pos, aa_wt, aa_mut):
        #Funcão que insere a mutação na sequência wildtype, dando a posição, o aa da sequência wildtype e o aa da sequência mutada.
        
        for i, val in enumerate(pos): 
            mut_seq = ""
            if self.sequence[int(val)-1] == aa_wt[i]: #primeiro elemento do tuplo (aa. wt, aa. mut)
                mut_seq += self.sequence[:int(val)-1]
                mut_seq += str(aa_mut[i]) #aminoacido de substituição
                mut_seq += self.sequence[int(val):]
                self.mut_sequence.append(str(mut_seq))
        return self.mut_sequence
        
#Criar péptidos:
    
    def create_peptide(self, protein_name, pos):
        #Função que cria a sequência peptídica para ser fragmentada e analisada pelo NetMHCpan, apartir do tamanho de péptido e do aa que sofre mutação.
     
        for tam in range(8,21):
            n = 0
            file = open(str(tam) + "_" + protein_name + ".fasta", "w")
            save_val = []
            for i, val in enumerate(pos):
                n += 1 
                if val not in save_val:
                    save_val.append(val)
                    wt = (str(self.sequence[val - tam : val + tam - 1]))
                    file.write(">WT" +str(n) + "\n" + wt + "\n \n")
                seq = self.mut_sequence[i]                  
                mut = (str(seq[val - tam : val + tam - 1]))
                file.write(">MUT" + str(n) + "\n" + mut + "\n \n")
        file.close()
        
        
if __name__=="__main__":
    
    def test1(): #IDH1
        s = sequence("49168486","IDH1_sequence.fasta")
        #s.protein_ncbi()
        s.protein_file()
        #s.protein_fasta("IDH1_sequence")
        s.induce_mutation([132],["R"],["H"]) 
        s.create_peptide("IDH1",[132])

    
    def test2(): #IDH2
        s = sequence("20141568","IDH2_sequence.fasta")
        #s.protein_ncbi()
        s.protein_file()
        #s.protein_fasta("IDH2_sequence")       
        s.induce_mutation([172,172,172],["R","R","R"],["G","K","M"]) 
        s.create_peptide("IDH2",[172,172,172])

    
    def test3(): #TP53
        s = sequence("823670074","TP53_sequence.fasta") 
        #s.protein_ncbi()
        s.protein_file()
        #s.protein_fasta("TP53_sequence")
        s.induce_mutation([143,175,248,273,281],["V","R","R","R","D"],["A","H","W","H","G"]) 
        s.create_peptide("TP53",[143,175,248,273,281])
    
    def test4(): #Histone H3-3
        s = sequence("55977062","H3_sequence.fasta") 
        #s.protein_ncbi()
        s.protein_file()
        #s.protein_fasta("H3_sequence")
        s.induce_mutation([28,35,35],["K","G","G"],["M","R","V"]) 
        s.create_peptide("H3",[28,35,35]) 

    test1()
    test2()
    test3()
    test4()
