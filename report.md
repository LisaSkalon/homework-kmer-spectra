## Отчет

Класс "Kmer_spectrum" позволяет рассчитывать спектр камеров заданной длины и заданного качества, визуализировать спектр и подсчитывать 
прикидочный размер генома как до отсечки шумовой части спектра, так и после.

Класс содержит методы:  
```python
analysis(self, in_file, k, q) # По входным данным (fastq файл, длина камера и качество) формирует словарь камер:частота встречаемости
```
```python
build(self) # Превращает словарь в np.array частота встречаемости - частота возникновения
```
```python
vizualize(self, xlim, ylim) # Из словаря и входных масштабов рисует график распределения камеров
```
```python
genome_size(self, z) # Рассчитывает примерный размер генома, учитывая шум и нет
```

Для файла ***test_kmer.fastq*** с размером камера 15 и отсечкой по качеству 35 был получен спектр:  

![plot](https://github.com/LisaSkalon/homework-kmer-spectra/blob/master/plot1.png)

Выбранная отсечка по шуму - 25.

Размер генома с шумом: 8321680.33333  
Размер без шума: 6805546.0  

