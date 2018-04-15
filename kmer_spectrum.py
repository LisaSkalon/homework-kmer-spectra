from Bio import SeqIO
from collections import defaultdict
import collections
import numpy as np
import matplotlib.pyplot as plt
import pylab
from matplotlib.ticker import NullFormatter, MultipleLocator
import argparse

# Парсинг
def parse_inputs():
     parser = argparse.ArgumentParser(description='Kmer spectra')
     parser.add_argument('-i', '--input' , help='Input file' , metavar='Str',
                    type=str, required=True)
     parser.add_argument('-k', '--kmersize', help = 'Kmer size', metavar = 'Int', type = int, default = 15)
     parser.add_argument('-q', '--quality', help = 'Nucleotide quality', metavar = 'Int', type = int, default = 35)
     parser.add_argument('-x', '--xlim', help = 'limit of x axis', metavar = 'Int', type = int, default = 150)
     parser.add_argument('-y', '--ylim', help = 'limit of y axis', metavar = 'Int', type = int, default = 150000)
     parser.add_argument('-z', '--granizza', help = 'limit of kmer spectra', metavar = 'Int', type = int, default = 25)
     
     args = parser.parse_args()
     return args.input, args.kmersize, args.quality, args.xlim, args.ylim, args.granizza
    
# Создаем класс с информацией о kmers.
    
class Kmer_spectrum:   
    
    def __init__(self):
        self.kmer_dict1 = defaultdict(int)
        self.kmer_list = []
        self.size_before = 0
        self.size_after = 0
        self.x = []
        self.y = []

    # Метод - анализирование спектра.
    def analysis(self, in_file, k, q):
        seq = ''
     
        with open(in_file, "r") as handle:
        # Читаем fasta файл, содержимое записываем в переменную.
            records = list(SeqIO.parse(handle, "fastq"))
        handle.close()
        
        # Проходимся по всем ридам и по каждому камеру, если у него хорошее 
        # качество, записываем его в словарь. Также в словарь записываются
        # частоты встречаемости каждого камера.
        for record in records:             
            seq = str(record.seq)
            seq_lng = len(seq)
            
            for index in range(seq_lng-k+1):
                
                qual = list(record.letter_annotations["phred_quality"][index:(index+k)])
                if all(item > q for item in qual):
                    current_kmer = seq[index:(index+k)]
                    self.kmer_dict1[current_kmer] +=1
    # Метод - построение спектра. Перевод словаря в массив.            
    def build(self):
        
        # Сортируем все значения встречаемости, считаем, сколько раз какое
        # значение встретилось с помощью каунтера.
        values = sorted(list(self.kmer_dict1.values()))
        counter=collections.Counter(values)
        # Создаем массив втречаемости и частоты.
        self.kmer_list = np.array([list(counter.values()), list(counter.keys())])
    # Метод - визуализайия спектра
    def vizualize(self, xlim, ylim):
        
        self.x = self.kmer_list[1]
        self.y = self.kmer_list[0]
       
        figure = pylab.figure(1)
        axes = figure.add_subplot (1, 1, 1)
        pylab.xlim (0, xlim)
        pylab.ylim (0, ylim)
        plt.plot(self.x,self.y, lw = 0.5)
        plt.scatter(self.x, self.y, lw = 1)
        
        minorLocator = MultipleLocator (base=5)
        axes.xaxis.set_minor_locator(minorLocator)
        axes.xaxis.set_minor_formatter(NullFormatter())
        
        plt.title( 'Kmer spectra', fontsize = 20,  loc='left')
        plt.xlabel('Occurrence')
        plt.ylabel(('Frequence'), fontsize = 12)
        
        plt.show()
        figure.savefig('plot1.pdf')
    # Метод - подсчет размера генома             
    def genome_size(self, z):

        summa = 0
        summa1 = 0
        # Подсчет без отсечки по шуму
        for i in range(len(self.kmer_list[0])):
            summa += (self.x[i]*self.y[i])
        self.size_before = summa/k
    
        # Подсчет с отсечкой по шуму
        f = self.x[z:]
        g = self.y[z:]
        
        for i in range(len(self.kmer_list[0][z:])):
            summa1 += (f[i]*g[i])
        self.size_after = summa1/k
        
        print('Genome size before cutting:', self.size_before)
        print('Genome size after cutting:', self.size_after)
 
# Получаем информацию по нашему классу. Выводит график и размер генома.       
if __name__ == '__main__':
    in_file, k, q, xlim, ylim, z = parse_inputs()
    
    spectra = Kmer_spectrum()
    spectra.analysis(in_file, k, q)
    spectra.build()
    spectra.vizualize(xlim, ylim)
    spectra.genome_size(z)
    
    
        
    

