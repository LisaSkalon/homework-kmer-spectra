# Third python homework

## Kmer spectra
This program draws kmer spectra and estimates genome size before and after noise cutting. 

## Usage and arguments
usage: kmer_spectrum.py [-h] -i Str [-k Int] [-q Int] [-x Int] [-y Int] [-z Int]  

required arguments:
 -i Str, --input Str   Input file  

optional arguments:  
  -h, --help            show this help message and exit    
  -k Int, --kmersize Int  
                        Kmer size  (def 15)  
  -q Int, --quality Int  
                        Nucleotide quality (def 35)   
  -x Int, --xlim Int    limit of x axis  (def 150)  
  -y Int, --ylim Int    limit of y axis  (def 150000)  
  -z Int, --granizza Int  
                        limit of kmer spectra  (def 25)  
                        
## Example 
input:  
python3 kmer_spectrum.py -i test_kmer.fastq

output:  

![plot](https://github.com/LisaSkalon/homework-kmer-spectra/blob/master/plot1.pdf)

Genome size before cutting: 8321680.33333  
Genome size after cutting: 6805546.0 





